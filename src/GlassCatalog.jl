module GlassCatalog

using YAML

const _dbFolder = joinpath(Pkg.dir("GlassCatalog"), "assets","refractiveindex.info-database","database")
const _libFile = joinpath(_dbFolder, "library.yml")

# Prepare for Julia 0.5 unified strings by overriding deprecated Base.String
typealias String UTF8String

immutable LibEntry
  name::String
  alias::String
  category::String
  catalias::String
  shelf::String
  path::String
end

string(e::LibEntry) = "{N:$name A:$alias C:$category CA:$catalias S:$shelf P:$path}"

#convert(::Type{Array{String,1}}, e::LibEntry) = [e.name, e.alias, e.category, e.catalias, e.shelf, e.path]

strip_html(s::AbstractString) = replace(s, r"<[^>]*>", "")

_library = Dict{String,LibEntry}()
_shelf_info = Dict{String,String}()
_paths = Dict{String,LibEntry}()

#= Form of the YAML library.yml file
 - SHELF: main
  name: "MAIN - simple inorganic materials"
  content:
    - DIVIDER: "Ag - Silver"
    - BOOK: Ag
      name: "Ag (Silver)"
      content:
        - PAGE: Rakic
          name: "Rakić 1998: n,k 0.2066-12.40 µm"
          path: "main/Ag/Rakic.yml"
        - PAGE: McPeak
          name: "McPeak et al. 2015: Thin film; n,k 0.3-1.7 µm"
          path: "main/Ag/McPeak.yml"
        ...
=#

# returns a LibEntry
function _parse_page(shelf, shelfdiv, book, bookdiv, bookname, page, pagename, path)
  shelf = lowercase(shelf)
  local spath = split(path, '/')
  local lc_shelfdiv = lowercase(shelfdiv)
  if shelf == "main"
    # shelfdiv form: "Al - Aluminium and aluminates"
    # parse as: "$div_formula - $div_name"
    m = match(r"^([^\-]+)\s\-\s(.+)$", shelfdiv)
    @assert m !== nothing && length(m.captures) >= 2 "Bad shelfdiv in $shelf: $shelfdiv"
    (div_formula, div_name) = m.captures

    # book form: "MgAl2O4"
    # bookname form: "MgAl<sub>2</sub>O<sub>4</sub> (Magnesium aluminate, spinel)"
    # parse as: "$bk_formula ($bk_name)"
    m = match(r"^(.+)\s\((.+)\)$", bookname)
    @assert m !== nothing && length(m.captures) >= 2 "Bad bookname in $shelf: $bookname"
    (bk_formula, bk_name) = m.captures

    # page form: "Tropf"
    # pagename form: "Tropf and Thomas 1991: n 0.35-5.5 µm"
    # leave unparsed

    # path form: "glass/schott/FK5HTi.yml"
    # parse as: ".../$(file_name).yml"
    @assert length(spath) > 0
    file_name = split(spath[end], '.')[1]

    entry_name = bookdiv == "" ? "$book-$file_name" : "$bookdiv-$file_name"

    # return a LibEntry(name, alias, category, catalias, shelf, path)
    return LibEntry(entry_name, "$bk_name - $pagename", div_formula, div_name, shelf, path)

  elseif shelf == "organic"
    # shelfdiv form: "Benzene and its derivatives"
    # leave unparsed

    # book form: "benzene"
    # bookname form: "C<sub>6</sub>H<sub>6</sub> (Benzene)"
    # parse as: "$bk_formula ($bk_name)"
    m = match(r"^(.+)\s\((.+)\)$", bookname)
    @assert m !== nothing && length(m.captures) >= 2 "Bad bookname in $shelf: $bookname"
    (bk_formula, bk_name) = m.captures

    # page form: "Moutzouris"
    # pagename form: "Moutzouris et al. 2013: n 0.450-1.551 µm"
    # leave unparsed

    # return a LibEntry(name, alias, category, catalias, shelf, path)
    return LibEntry("$book-$page", "$bk_name - $pagename", shelfdiv, shelfdiv, shelf, path)

  elseif shelf == "glass" && startswith(lc_shelfdiv, "optical")
    # shelfdiv form: "Optical glass - SCHOTT"
    # leave unparsed

    # book form: "SCHOTT-FK"
    # bookname form: "SCHOTT - FK (Fluorite crown)"
    #   OR "SCHOTT - Inquiry glasses"
    # parse as: "$bk_mfg - $bk_class ($bk_classname)"
    #   OR "$bk_mfg - $bk_classname"
    if ismatch(r"^(.+)\s-\s(.+)\s\((.+)\)$", bookname)
      m = match(r"^(.+)\s-\s(.+)\s\((.+)\)$", bookname)
      @assert m !== nothing && length(m.captures) >= 3 "Bad bookname in $shelf: $bookname"
      (bk_mfg, bk_class, bk_classname) = m.captures
    else
      m = match(r"^(.+)\s-\s(.+)$", bookname)
      @assert m !== nothing && length(m.captures) >= 2 "Bad bookname in $shelf: $bookname"
      (bk_mfg, bk_classname) = m.captures
      m = nothing #match(r"^(.+) glasses$", bk_classname)
      if m === nothing
        bk_class = bk_classname
      else
        (bk_class) = m.captures
      end
    end

    # page form: "FK5HTi"
    # pagename form: "FK5HTi"
    # leave unparsed

    # return a LibEntry(name, alias, category, catalias, shelf, path)
    return LibEntry("$page-$bk_mfg", "$pagename - $bk_mfg", "$bk_mfg - $bk_class", "$shelfdiv - $bk_classname", shelf, path)


  elseif shelf == "glass" && startswith(lc_shelfdiv, "popular")
    # popular glasses are listed twice. Modify the entry name to start with popular
    # return a LibEntry(name, alias, category, catalias, shelf, path)
    return LibEntry("popular-$book-$page", "$bookname - $pagename (popular)", shelfdiv, shelfdiv, shelf, path)

  elseif shelf == "glass"
          # shelfdiv form: "Multi purpose glass"
          # leave unparsed

          # book form: "soda-lime"
          # bookname form: "Soda lime glass"
          #   OR "SCHOTT - multiple purpose"
          # parse as: "$bk_mfg - $bk_class ($bk_classname)"
          #   OR "$bk_mfg - $bk_classname"
          if ismatch(r"^(.+)\s-\s(.+)$", bookname)
            m = match(r"^(.+)\s-\s(.+)$", bookname)
            @assert m !== nothing && length(m.captures) >= 2 "Bad bookname in $shelf: $bookname"
            (bk_mfg, bk_class) = m.captures
            bk_classname = bk_class
          else
            bk_mfg = ""
            bk_class = bookname
            bk_classname = bk_class
          end

          # page form: "Rubin-clear"
          #   OR "ZERODUR"
          # leave unparsed

          if length(bk_mfg) == 0
            # return a LibEntry(name, alias, category, catalias, shelf, path)
            return LibEntry("$book-$page", "$bookname - $pagename", shelfdiv, shelfdiv, shelf, path)
          else
            # return a LibEntry(name, alias, category, catalias, shelf, path)
            return LibEntry("$page-$bk_mfg", "$pagename - $bk_mfg", "$bk_mfg - $bk_class", "$shelfdiv - $bk_classname", shelf, path)
          end



  else
    if ismatch(r"^(.+)\s\((.+)\)$", bookname)
      m = match(r"^(.+)\s\((.+)\)$", bookname)
      @assert m !== nothing && length(m.captures) >= 2 "Bad bookname in $shelf: $bookname"
      (bk_formula, bk_name) = m.captures

      # return a LibEntry(name, alias, category, catalias, shelf, path)
      return LibEntry("$bk_formula-$page", "$bk_name - $pagename", shelfdiv, shelfdiv, shelf, path)
    else
      if shelfdiv==""
        shelfdiv = _shelf_info[shelf]
      end
      # return a LibEntry(name, alias, category, catalias, shelf, path)
      return LibEntry("$book-$page", "$bookname - $pagename", shelfdiv, shelfdiv, shelf, path)
    end
  end
end

function _parse_library()
  lib = YAML.load_file(_libFile)
  entry_count::Int = 0

  for shelf_d in lib
    (shelf, shelfname) = (shelf_d["SHELF"], shelf_d["name"])
    _shelf_info[shelf] = shelfname
    #if lowercase(shelf) != "main"
    #  continue
    #end
    shelfdiv::String = ""
    for book_d in shelf_d["content"]
      if haskey(book_d,"DIVIDER")
        shelfdiv = book_d["DIVIDER"]
      else
        (book, bookname) = (book_d["BOOK"], strip_html(book_d["name"]))
        bookdiv::String = ""
        for page_d in book_d["content"]
          if haskey(page_d,"DIVIDER")
            bookdiv = page_d["DIVIDER"]
          end
          if haskey(page_d, "PAGE")
            (page, pagename, path) = (page_d["PAGE"], strip_html(page_d["name"]), page_d["path"])
            entry = _parse_page(shelf, shelfdiv, book, bookdiv, bookname, page, pagename, path)
            if haskey(_library, entry.name)
              oldentry = _library[entry.name]
              println("WARNING: $(entry.name) in $(entry.catalias) already in library under $(oldentry.catalias)\nOLD=$oldentry\nNEW=$entry")
            elseif haskey(_paths, path)
                oldentry = _paths[path]
                #println("WARNING: $path for $(entry.name) already in library under $(oldentry.name)\nOLD=$oldentry\nNEW=$entry")
            else
              _library[entry.name] = entry
              _paths[path] = entry
              entry_count += 1
              if entry_count > 10000
                return
              end
            end
          end
        end
      end
    end
  end
end


using DataFrames

function as_dataframe()
  _parse_library()

  a_key = String[]
  a_name = String[]
  a_alias = String[]
  a_category = String[]
  a_catalias = String[]
  a_shelf = String[]
  a_path = String[]
  for k in sort(collect(keys(_library)))
    e = _library[k]
    push!(a_key, k)
    push!(a_name, e.name)
    push!(a_alias, e.alias)
    push!(a_category, e.category)
    push!(a_catalias, e.catalias)
    push!(a_shelf, e.shelf)
    push!(a_path, e.path)
  end
  return DataFrame(key=a_key, name=a_name, alias=a_alias, catagory=a_category, catalias=a_catalias, shelf=a_shelf, path=a_path)
end

writetable(joinpath(homedir(),"Devel/sandbox/julia/df.csv"), as_dataframe())

end # module GlassCatalog
