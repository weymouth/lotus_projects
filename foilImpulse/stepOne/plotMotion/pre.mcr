#!MC 1120
# Created by Tecplot 360 build 11.3.29.563
$!READDATASET  ' "../fort.300" "../fort.301" "../fort.302" "../fort.303" "../fort.304" "../fort.305" "../fort.306" "../fort.307" "../fort.308" "../fort.309" "../fort.310" "../fort.311" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "p"'
$!GLOBALCONTOUR 1  VAR = 4

$!WRITEDATASET  "./pre.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INCLUDEDATASHARELINKAGE = YES
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [1-|NUMZONES|]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9

$!VarSet |first_line| = |NUMZONES|
$!VarSet |first_line| += 1

$!SLICELAYERS SHOW = YES
$!SLICEATTRIBUTES 1  SHOWGROUP = NO
$!SLICEATTRIBUTES 3  SHOWGROUP = YES
$!SLICEATTRIBUTES 3  PRIMARYPOSITION{Z = 0.71428}

$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'Extract Over Time'
  COMMAND = 'ExtractSliceOverTime'

$!WRITEDATASET  "./pre_Z0p714.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INCLUDEDATASHARELINKAGE = YES
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [|first_line|-|NUMZONES|]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9

$!VarSet |first_line| = |NUMZONES|
$!VarSet |first_line| += 1

$!ISOSURFACELAYERS SHOW = YES
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 0.3

$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'Extract Over Time'
  COMMAND = 'ExtractIsoSurfaceOverTime'

$!WRITEDATASET  "./pre_P0p3.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [|first_line|-|NUMZONES|]
  VARPOSITIONLIST =  [1-3]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9

$!VarSet |first_line| = |NUMZONES|
$!VarSet |first_line| += 1

$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = -0.3

$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'Extract Over Time'
  COMMAND = 'ExtractIsoSurfaceOverTime'

$!WRITEDATASET  "./pre_Pm0p3.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [|first_line|-|NUMZONES|]
  VARPOSITIONLIST =  [1-3]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9
