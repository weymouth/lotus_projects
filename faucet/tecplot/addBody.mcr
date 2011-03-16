#!MC 1120
# Created by Tecplot 360 build 11.3.29.563

$!Varset |blocks| = 12
$!Varset |body| = 100
$!Varset |press| = 300

$!NEWLAYOUT 

$!Loop |blocks|

$!Varset |bnum| = (|Loop|+|body|-1)
$!Varset |pnum| = (|Loop|+|press|-1)

$!READDATASET  ' "fort.|bnum|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "p" "dis"'

$!Varset |bzone| = |NUMZONES|

$!READDATASET  ' "fort.|pnum|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "p" "dis"'

$!Varset |steps| = (|NUMZONES|-|bzone|)

$!Loop |steps|
$!Varset |pzone| = (|bzone|+|Loop|)
$!ALTERDATA  [|pzone|]
  EQUATION = '{dis} = {p}[|bzone|]'
$!Endloop

$!Varset |first| = (|bzone|+1)
$!WRITEDATASET  "./addBody.|pnum|.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [|first|-|NUMZONES|]
  INCLUDEDATASHARELINKAGE = YES
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9

$!Endloop