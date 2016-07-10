using YAML

dbfolder = Pkg.dir("GlassCatalog","assets","refractiveindex.info-database","database")
dbfile = joinpath(dbfolder, "library.yml")
library = YAML.load_file(dbfile)

using DataArrays, DataFrames

shelves = Dict{UTF8String,UTF8String}()
books = Dict{UTF8String,UTF8String}()
pages = Dict{UTF8String,UTF8String}()

libraryData = DataFrame(shelf=UTF8String[], subset=UTF8String[], book=UTF8String[], page=UTF8String[], path=UTF8String[])

for shelf in library
  shelfName::UTF8String = shelf["SHELF"]
  shelfDesc::UTF8String = shelf["name"]
  setindex!(shelves, shelfName, shelfDesc)
  subset::UTF8String = ""
  for book in shelf["content"]
    if haskey(book,"DIVIDER")
      subset = book["DIVIDER"]
    else
      bookName::UTF8String = book["BOOK"];
      bookDesc::UTF8String = book["name"]
      setindex!(books, bookName, bookDesc)
      for page in book["content"]
        if haskey(page, "PAGE")
          pageName::UTF8String = page["PAGE"]
          pageDesc::UTF8String = page["name"]
          pagePath::UTF8String = page["path"]
          setindex!(pages, pageName, pageDesc)
          df = DataFrame(shelf=[shelfName], subset=[subset], book=[bookName], page=[pageName], path=[pagePath])
          append!(libraryData, df)
        end
      end
    end
  end
end
