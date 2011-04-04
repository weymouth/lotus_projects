#!MC 1120
# Created by Tecplot 360 build 11.3.29.563

$!Varset |blocks| = 12
$!Varset |velo| = 200
$!Varset |lam2| = 400

$!Loop |blocks|

$!NEWLAYOUT 

$!Varset |vnum| = (|Loop|+|velo|-1)
$!Varset |lnum| = (|Loop|+|lam2|-1)

$!READDATASET  ' "fort.|vnum|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "u" "v" "w" "vort"'

$!Varset |steps| = |NUMZONES|

$!READDATASET  ' "fort.|lnum|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "u" "v" "w" "vort" "p" '

$!Loop |steps|
$!Varset |lzone| = (|steps|+|Loop|)
$!ALTERDATA  [|lzone|]
  EQUATION = '{vort} = ddx({v}[|Loop|])-ddy({u}[|Loop|])'
$!Endloop

$!RENAMEDATASETVAR 
  VAR = 8
  NAME = 'lam2'

$!ALTERDATA 
  EQUATION = '{dis} = 1'
$!ALTERDATA  [1]
  EQUATION = '{dis} = 0.15-{x}*sin(pi/18)+{y}*cos(pi/18)'
$!ALTERDATA  [1]
  EQUATION = '{dis} = max({dis},0.4-{z})	'
$!ALTERDATA 
  EQUATION = '{dis} = {dis}[1]'

$!Varset |first| = (|steps|+1)
$!WRITEDATASET  "./Lam2AndVort.|lnum|.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  ASSOCIATELAYOUTWITHDATAFILE = NO
  INCLUDEDATASHARELINKAGE = YES
  ZONELIST =  [|first|-|NUMZONES|]
  VARPOSITIONLIST =  [1-3,7-9]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9

$!Endloop