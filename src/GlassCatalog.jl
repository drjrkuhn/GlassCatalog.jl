module GlassCatalog

using YAML

const _dbFolder = joinpath(Pkg.dir("GlassCatalog"), "assets","refractiveindex.info-database","database")
const _libFile = joinpath(_dbFolder, "library.yml")

# Prepare for Julia 0.5 unified strings by overriding deprecated Base.String
typealias String UTF8String

immutable LibEntry
  path::String    # YAML file for the material (page)
  description::String # description of material
  subgroup::String  # subset of materials within group
  group::String   # group of materials (book)
  category::String    # overall group category (shelf)
end

_library = Dict{String,LibEntry}()
_groupInfo = Dict{String,String}()
_categoryInfo = Dict{String,String}()

function _parse_library()
  lib = YAML.load_file(_libFile)
  count::Int = 0

  for shelf in lib
    shelfName::String = shelf["SHELF"]
    shelfDesc::String = shelf["name"]
    _categoryInfo[shelfName] = shelfDesc
    subset::String = ""
    for book in shelf["content"]
      if haskey(book,"DIVIDER")
        subset = book["DIVIDER"]
      else
        bookName::String = book["BOOK"];
        bookDesc::String = book["name"]
        _groupInfo[bookName] = bookDesc
        for page in book["content"]
          if haskey(page, "PAGE")
            pageName::String = page["PAGE"]
            pageDesc::String = page["name"]
            pagePath::String = page["path"]
            entry = LibEntry(pagePath, pageDesc, subset, bookName, shelfName)
            if haskey(_library, pageName)
              error("$pageName already in library")
            end
            println(">> $pageName: $entry")
            _library[pageName] = entry
            count+=1
            if count > 20
              return
            end
          end
        end
      end
    end
  end
end

end # module GlassCatalog
