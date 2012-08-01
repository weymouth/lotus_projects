#!MC 1120
# Created by Tecplot 360 build 11.3.29.563
$!Varset |blocks| = 12
$!Varset |lamb| = 400
$!Varset |vort| = 200

$!NEWLAYOUT 

$!VarSet |last| = 0

$!Loop |blocks|

$!Varset |vnum| = (|Loop|+|vort|-1)
$!Varset |lnum| = (|Loop|+|lamb|-1)

$!READDATASET  ' "fort.|vnum|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "u" "v" "w" "p"'

$!Varset |steps| = (|NUMZONES|-|last|)

$!READDATASET  ' "fort.|lnum|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "u" "v" "w" "p"'

$!Loop |steps|
$!Varset |vzone| = (|last|+|Loop|)
$!Varset |lzone| = (|vzone|+|steps|)
$!ALTERDATA  [|vzone|]
  EQUATION = '{p} = {p}[|lzone|]'
$!Endloop

$!Varset |first| = (|last|+|steps|+1)
$!DELETEZONES [|first|-|NUMZONES|]
$!VarSet |last| = |NUMZONES|

$!Endloop

$!GLOBALCONTOUR 1  VAR = 7
$!ISOSURFACELAYERS SHOW = YES
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = -10.0
$!ISOSURFACEATTRIBUTES 1  OBEYSOURCEZONEBLANKING = YES

$!RUNMACROFUNCTION  "IJKBlank"

$!VarSet |first_zone| = (|NUMZONES|+1)
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'Extract Over Time'
  COMMAND = 'ExtractIsoSurfaceOverTime'

#$!DELETEVARS [4,5]
#$!CREATEMIRRORZONES 
#  SOURCEZONES =  [|first_zone|-|NUMZONES|]
#  MIRRORVAR = 'Y'

$!WRITEDATASET  "lam2.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INCLUDEDATASHARELINKAGE = YES
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [|first_zone|-|NUMZONES|]
  VARPOSITIONLIST =  [1-6]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9
