using YAML, DataFrames

const _dbFolder = joinpath(Pkg.dir("GlassCatalog"), "assets","refractiveindex.info-database","database")
const _libFile = joinpath(_dbFolder, "library.yml")

# Prepare for Julia 0.5 unified strings by overriding deprecated Base.String
typealias String UTF8String

strip_html(s::AbstractString) = replace(s, r"<[^>]*>", "")

_library = DataFrame(Name=String[], Alias=String[], Category=String[], CatAlias=String[], Shelf=String[], Path=String[])
_shelf_info = Dict{String,String}()

#= Form of the YAML library.yml file
 - SHELF: main
  name: "MAIN - simple inorganic materials"
  content:
    - DIVIDER: "Ag - Silver"
    - BOOK: Ag
      name: "Ag (Silver)"
      content:
        - DIVIDER: "optional sub-grouping of pages"
        - PAGE: Rakic
          name: "Rakić 1998: n,k 0.2066-12.40 µm"
          path: "main/Ag/Rakic.yml"
        - PAGE: McPeak
          name: "McPeak et al. 2015: Thin film; n,k 0.3-1.7 µm"
          path: "main/Ag/McPeak.yml"
        ...
=#

# returns an array of strings representing a row in the dataframe
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

    entry_name = isempty(bookdiv) ? "$book-$file_name" : "$bookdiv-$file_name"

    # return (name, alias, category, catalias, shelf, path)
    return (entry_name, "$bk_name - $pagename", div_formula, div_name, shelf, path)

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

    # return (name, alias, category, catalias, shelf, path)
    return ("$book-$page", "$bk_name - $pagename", shelfdiv, shelfdiv, shelf, path)

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

    # return (name, alias, category, catalias, shelf, path)
    return ("$page-$bk_mfg", "$pagename - $bk_mfg", "$bk_mfg - $bk_class", "$shelfdiv - $bk_classname", shelf, path)


  elseif shelf == "glass" && startswith(lc_shelfdiv, "popular")
    # popular glasses are listed twice. Modify the entry name to start with popular
    # return (name, alias, category, catalias, shelf, path)
    return ("popular-$book-$page", "$bookname - $pagename (popular)", shelfdiv, shelfdiv, shelf, path)

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
            # return (name, alias, category, catalias, shelf, path)
            return ("$book-$page", "$bookname - $pagename", shelfdiv, shelfdiv, shelf, path)
          else
            # return (name, alias, category, catalias, shelf, path)
            return ("$page-$bk_mfg", "$pagename - $bk_mfg", "$bk_mfg - $bk_class", "$shelfdiv - $bk_classname", shelf, path)
          end



  else
    if ismatch(r"^(.+)\s\((.+)\)$", bookname)
      m = match(r"^(.+)\s\((.+)\)$", bookname)
      @assert m !== nothing && length(m.captures) >= 2 "Bad bookname in $shelf: $bookname"
      (bk_formula, bk_name) = m.captures

      # return (name, alias, category, catalias, shelf, path)
      return ("$bk_formula-$page", "$bk_name - $pagename", shelfdiv, shelfdiv, shelf, path)
    else
      if isempty(shelfdiv)
        shelfdiv = _shelf_info[shelf]
      end
      # return (name, alias, category, catalias, shelf, path)
      return ("$book-$page", "$bookname - $pagename", shelfdiv, shelfdiv, shelf, path)
    end
  end
end

function _parse_library()
  # load the .yml file as a Dict hierarchy
  lib_d = YAML.load_file(_libFile)

  for shelf_d in lib_d
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
            push!(_library, entry)
          end
        end
      end
    end
  end
end

using DataFrames

function as_dataframe()
  _parse_library()
  return _library
end

writetable(joinpath(homedir(),"Devel/sandbox/julia/df.csv"), as_dataframe())
