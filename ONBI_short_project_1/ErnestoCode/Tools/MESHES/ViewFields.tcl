package require vtk
package require vtkinteraction

proc ReadAttributes { } {
  .att.m delete 0 end
  .att.m add command -label "None" -command "ShowNone"
  for { set i 0 } { $i < [ [ OBJETO GetPointData] GetNumberOfArrays ] } { incr i } {
    set att_name [[[ OBJETO GetPointData] GetArray  $i ] GetName ]
    if { $att_name != "__actual__" } {
      set att_nofc [[[ OBJETO GetPointData] GetArray  $i ] GetNumberOfComponents ]

      if        { $att_nofc == 1 } {
        if { [lsearch  [image names] "NodeScalar"] < 0 } {
          image create photo NodeScalar -data R0lGODlhEAAQAPf/AP///+YLJ+UJHuELJuALJtwKJtcLJdYOLtYNLdQNLdMNLdMMLcsQNcsPNckIG8QSPcQRPcQRPMMSPcMRPbwURbwURLwTRbwTRLUWTLUVTLQWTLQVTLMZV7IYV7EYV68cYq8YVq0YU60XVK0XU6wYVKwXVKobXqkaXqcaXKYeaaUaW6UZXJ8cZJ4cZJ4cY50dY50cY5wfb5wcYZgfbJgXVJYea5YeapUfa5Uea5MeaZIidpIhdpIeaY8hc44hco4gcoogb4kkfIgje4cieoYjeoYieoYYVIQYU4MmhoImhoAlgn8lgn8kgX8TRX4kgn0pjnwkfnsFEXknincniXYrmHYnh3Umh3Eqk3AuonAqknApkW8pkW8dZW8cZW4pjmstnWosnWksmWgsmWgsmGgrmGcrmWQvpmQvpWMupGIuo2IPNWAuoGAtoF8ojF8TQl4je10je10EDVkwqVkwqFkvqFgwqFIPMlEysFEyr1Exr1ESQE4hdEk0t0kYVkgVSUUCB0I2v0I2vkE2v0E2vjwaWzs6zDo5yDo5xjo4xjYUSDM81DM70jI80TI7zi00uCg5yScKIyU0uCAHGR4HFx4GFBwIHRoJHhkHGw8TRA0ROwwDCAsEDwoAAAUBBAQBBAMBAgIEDgIAAwEDCQEAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACwAAAAAEAAQAAAIywA7deLkIMCBBg8oYACBgsUkgQLjGGSQEIOIhjH0CPwTxcDBihdZ1OgxRE8lARMrcviQYuQQJFe8FPiokEMTO26ABFlS5coXmhg4qIEICYqUK1vMIFQIggYniJ36UElqBqQJI1A7JTJjBg8epg1rPIS4x4wcPHwWhu3BpVInUITMogUU0iWSKm/ayOUD6JAJkSSXTPG5F1AhRUZcCt7ypYxXvoUaRep0iQuSKVvKdJ1byFEmUAIvYS6jee6hRpiyviGNlu8hRY9AdwoIADs=
        }
        .att.m add command -label $att_name          -command "SetAttributesOnNodes $att_name 0"  -image NodeScalar  -compound left
      } elseif  { $att_nofc == 3 } {
        if { [lsearch  [image names] "NodeVector"] < 0 } {
          image create photo NodeVector -data R0lGODlhEAAQANX/AJCQkPf/9/b29unp6eTk5N7e3srKysXFxcP/w7T/tLHdsbGxsa+vr5zAnJv/m5CQkIb/hoGBgX//f3l5eXb/dm5ubm3/bWP/Y2JiYmJWYlxcXFRUVE5OTkn/ST//Pz8/Pyz/LCD/IBsbGxn/GQoKCgn/CQD/AAD2AADfAAC/AABjAABXAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAAAZVQIBwSCwaj8gjIaJJCgsTEaviXLCupIJz8LlunIAChsNiOAsbA+AhSKIPSYSFosjAkYmSCbVqOD0mJyopF04ggSgmHU4SJo4mDk4BFCMhEGBCAQFgQQA7
        }
        .att.m add command -label $att_name          -command "DrawVectorsOnNodes $att_name"      -image NodeVector  -compound left

        if { [lsearch  [image names] "NodeScalarC"] < 0 } {
          image create photo NodeScalarC -data R0lGODlhEAAQAID/AMDAwAAAACH5BAEAAAAALAAAAAAQABAAAAIqBIKGa735QlPvVBaprFZO7X0fJ45J2aFTpnbnerIOx5g17MngHPI9+igAADs=
        }
        .att.m add command -label "    $att_name (x)" -command "SetAttributesOnNodes $att_name 0"  -image NodeScalarC -compound right
        .att.m add command -label "    $att_name (y)" -command "SetAttributesOnNodes $att_name 1"  -image NodeScalarC -compound right
        .att.m add command -label "    $att_name (z)" -command "SetAttributesOnNodes $att_name 2"  -image NodeScalarC -compound right
      } elseif  { $att_nofc == 9 } {
        if { [lsearch  [image names] "NodeTensor"] < 0 } {
          image create photo NodeTensor -data R0lGODlhEAAQAMT/AMDAwP+rq/+IiP9/f/9GRv89Pf8yMv8qKv8dHf8UFP8AAPIBDdgCJ77B/7y//6yw/56j/3d9/251/2Vs/zxF/zoY1Co0/yQK2xsm/xEc/w0L8gEN/wAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAAAVPICCOZGmOQoIo53goMFy0CmfHdH0rxgkzNo5iYArAFpeNstIwEWIapVRSSsSkWAspQJIGlS2H8rt5tDAb8qYFmGCVFPYj82aLIhs0xG4KAQA7EAAQAAAGZ0CAcEgMEI/DACICmQyQw8dkOm1ABxOXllpAKqZa19ThnWbCkwTyML2YWHASB8mgquD4EFEQoeL/JUQHAwoFGnhhcFAAHXCJLB6LJyyPLBSLIn9wI4sAHimanQAUICyTH6JC
        }
        .att.m add command -label $att_name -command "DrawTensorsOnNodes $att_name"               -image NodeTensor  -compound left

        if { [lsearch  [image names] "NodeScalarC"] < 0 } {
          image create photo NodeScalarC -data R0lGODlhEAAQAID/AMDAwAAAACH5BAEAAAAALAAAAAAQABAAAAIqBIKGa735QlPvVBaprFZO7X0fJ45J2aFTpnbnerIOx5g17MngHPI9+igAADs=
        }
        for { set c 0 } { $c < 9 } {incr c } {
          .att.m add command -label "    $att_name ([expr $c+1])" -command "SetAttributesOnNodes $att_name $c"   -image NodeScalarC -compound right
        }
      } else {
        if { [lsearch  [image names] "NodeScalarC"] < 0 } {
          image create photo NodeScalarC -data R0lGODlhEAAQAID/AMDAwAAAACH5BAEAAAAALAAAAAAQABAAAAIqBIKGa735QlPvVBaprFZO7X0fJ45J2aFTpnbnerIOx5g17MngHPI9+igAADs=
        }
        for { set c 0 } { $c < $att_nofc } {incr c } {
          .att.m add command -label "$att_name ([expr $c+1])" -command "SetAttributesOnNodes $att_name $c"      -image NodeScalarC -compound left
        }
      }
    }
  }

  for { set i 0 } { $i < [ [ OBJETO GetCellData] GetNumberOfArrays ] } { incr i } {
    set att_name [[[ OBJETO GetCellData] GetArray  $i ] GetName ]
    if { $att_name != "__actual__" } {
      set att_nofc [[[ OBJETO GetCellData] GetArray  $i ] GetNumberOfComponents ]

      if        { $att_nofc == 1 } {
        if { [lsearch  [image names] "CellScalar"] < 0 } {
          image create photo CellScalar -data R0lGODlhEAAQANX/AP8AAPQAAO8F3O0AANsAANMAAMgAAMAAAL8AAKwAAKEAAEYAADsAAC4AACcAABoHBxUHBg8AAQoBCgcHBgAa/wAX+wAV9wAQpAAOlAANoQAMfwAKZQAJagAENQADLwABYwABWAAADwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAIALAAAAAAQABAAAAZ3QIFwSCwKI4+kcumQRISMgXQ6DRgaE+ECwO12CQ0OSOv1FhwbyodcBpzTlLFgWz44NJRKRU7vHh54eWpsXAgQGhWCe4QIEYiKcWwJERgUkBVrcwAJEpV6kHwKEheWiaZ5mREhGXoWFa6wFnITIR0et7i5IUa8REEAOwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACwAAAAAEAAQAAAHl4ADCAyEhYYMCAgDGjUwjo+QLTg1GiUqBJiZmQIOKzoaIyYBo6SkBytLUR+hpaUJLUo7UB8loq0Br7E7qrW3Dy1IO0BAvLakDzDBwrKrxgEQMUhAy8S0xhA10tS7Hr0BFzVEO9tAs70XOOLD26qhGThD49PzwrOMVUbDQkD7/ULtOqo0cUKwoMEqGgpQoFChocOHFSQUCAQAOw==
        }
        .att.m add command -label $att_name          -command "SetAttributesOnCells $att_name 0"  -image CellScalar  -compound left
      } elseif  { $att_nofc == 3 } {
        if { [lsearch  [image names] "CellVector"] < 0 } {
          image create photo CellVector -data R0lGODlhEAAQAMT/AMDAwPDw8O3t7eHh4dzc3NXV1c/Pz7a2tqqqqqCgoJWVlY+Pj4eHh39/f3V1dWxsbGBgYFtbW1RUVE9PT0RERD09PTMzMywsLCUlJRcXFw4ODgUFBQAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAQAViICCOZGkCQbRt4zORK0semgaNBhadgCONMR7BcFoBGBtF70caXCwDkezBJDQaDseDKpFMJhEEbywNlgyZFxBgFCU2Cx9vpREvTYIKxtCWjwgYFAEmfgUrE10SW1kOF2pkJSEAOwyVy6LiVSJCi41YjGMFRAGQmeQCOi4nswiAjFTGHEK2UdD5aqrFAYqliMqHBCopBEIzQhplQwYsJwJGfwosHUhyCEolIJkaGhedKngGFBSdm5kgJSUdcEhIQQA7
        }
        .att.m add command -label $att_name          -command "DrawVectorsOnCells $att_name"      -image CellVector  -compound left

        if { [lsearch  [image names] "CellScalarC"] < 0 } {
          image create photo CellScalarC -data  R0lGODlhEAAQAMT/AMDAwMDAwNsAANMAAMgAAMAAAL8AAKwAAKEAAEYAADsAAC4AACcAABoHBxUHBg8AAQoBCgcHBgAQpAAOlAANoQAMfwAKZQAJagAENQADLwABYwABWAAADwAAAAAAAAAAACH5BAEAAAEALAAAAAAQABAAQAVdYCCOZCkmQAoYTqWm2/m+A2MB2olA0vzmAYXvRVhEZLMCwwVDHh6TITDCwWSu2CxH9Gh4v2AG5IEcChaXWADlq90Aaraq0GCmgHLDw66KAw4QUUNqDxwUQ30miiQhADs7
        }
        .att.m add command -label "    $att_name (x)" -command "SetAttributesOnCells $att_name 0"   -image CellScalarC  -compound right
        .att.m add command -label "    $att_name (y)" -command "SetAttributesOnCells $att_name 1"   -image CellScalarC  -compound right
        .att.m add command -label "    $att_name (z)" -command "SetAttributesOnCells $att_name 2"   -image CellScalarC  -compound right
      } elseif  { $att_nofc == 9 } {
        if { [lsearch  [image names] "CellTensor"] < 0 } {
          image create photo CellTensor -data  R0lGODlhEAAQAMT/AMDAwP+9vf+0tP+rq/+np/+Wlv+Ojv+Cgv9DQ/86Ov81Nf8dHf8UFP8NDf8CAtUAAGUAADIAAC0AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAQAVLICCOZClOgKOuLDKiwAQLDeMsTnECxhAYJdlLaDLBDKykgwCQKJSORyQoq1qtJ1gxGNNuYy8I1AUWQFcHkeSsYpAmh2Qi+/7u7MUQADsAAAAAEAAQAAAGW0CAcEgsCkUioxKATC6JSWTRADAICsejU6HBbC6bRvbo2ZjPG8qAKXSg3xsGWwFHR46G+pkyBOnNGEQgFnUdH0UiEW8VVFJDIggOCwkOa49ObE+ZTE2dnpiaREEAOw==
        }
        .att.m add command -label $att_name -command "DrawTensorsOnCells $att_name"               -image CellTensor  -compound left

        if { [lsearch  [image names] "CellScalarC"] < 0 } {
          image create photo CellScalarC -data  R0lGODlhEAAQAMT/AMDAwMDAwNsAANMAAMgAAMAAAL8AAKwAAKEAAEYAADsAAC4AACcAABoHBxUHBg8AAQoBCgcHBgAQpAAOlAANoQAMfwAKZQAJagAENQADLwABYwABWAAADwAAAAAAAAAAACH5BAEAAAEALAAAAAAQABAAQAVdYCCOZCkmQAoYTqWm2/m+A2MB2olA0vzmAYXvRVhEZLMCwwVDHh6TITDCwWSu2CxH9Gh4v2AG5IEcChaXWADlq90Aaraq0GCmgHLDw66KAw4QUUNqDxwUQ30miiQhADs7
        }
        for { set c 0 } { $c < 9 } {incr c } {
          .att.m add command -label "    $att_name ([expr $c+1])" -command "SetAttributesOnCells $att_name $c"   -image CellScalarC  -compound right
        }
      } else {
        if { [lsearch  [image names] "CellScalarC"] < 0 } {
          image create photo CellScalarC -data  R0lGODlhEAAQAMT/AMDAwMDAwNsAANMAAMgAAMAAAL8AAKwAAKEAAEYAADsAAC4AACcAABoHBxUHBg8AAQoBCgcHBgAQpAAOlAANoQAMfwAKZQAJagAENQADLwABYwABWAAADwAAAAAAAAAAACH5BAEAAAEALAAAAAAQABAAQAVdYCCOZCkmQAoYTqWm2/m+A2MB2olA0vzmAYXvRVhEZLMCwwVDHh6TITDCwWSu2CxH9Gh4v2AG5IEcChaXWADlq90Aaraq0GCmgHLDw66KAw4QUUNqDxwUQ30miiQhADs7
        }
        for { set c 0 } { $c < $att_nofc } {incr c } {
          .att.m add command -label "$att_name ([expr $c+1])" -command "SetAttributesOnCells $att_name $c"      -image CellScalarC  -compound left
        }
      }
    }
  }
}

proc ShowNone { } {
  MAPPER ScalarVisibilityOff
  if { [info command VECTORSACTOR] != "" } { VECTORSACTOR VisibilityOff  }
  if { [info command TENSORSACTOR] != "" } { TENSORSACTOR VisibilityOff  }
  VENTANA Render
}

proc SetAttributesOnCells { att { c 0 } } {
  if { [info command VECTORSACTOR] != "" } { VECTORSACTOR VisibilityOff  }
  if { [info command TENSORSACTOR] != "" } { TENSORSACTOR VisibilityOff  }

  if { $c == 0 } {
    [OBJETO GetCellData] SetActiveScalars   $att
  } else {
    vtkDoubleArray ATT
      ATT SetName __actual__
      ATT SetNumberOfComponents 1
      ATT SetNumberOfTuples [OBJETO GetNumberOfCells]
    ATT CopyComponent 0 [ [OBJETO GetCellData] GetArray $att] $c

    [OBJETO GetCellData] AddArray           ATT
    [OBJETO GetCellData] SetActiveScalars   __actual__
    ATT Delete
  }

  MAPPER SetScalarModeToUseCellData
  MAPPER ScalarVisibilityOn
  eval MAPPER SetScalarRange [[ [OBJETO GetCellData] GetScalars ] GetRange  ]
  VENTANA Render
}

proc SetAttributesOnNodes { att { c 0 } } {
  if { [info command VECTORSACTOR] != "" } { VECTORSACTOR VisibilityOff  }
  if { [info command TENSORSACTOR] != "" } { TENSORSACTOR VisibilityOff  }

  if { $c == 0 } {
    [OBJETO GetPointData] SetActiveScalars   $att
  } else {
    vtkDoubleArray ATT
      ATT SetName __actual__
      ATT SetNumberOfComponents 1
      ATT SetNumberOfTuples [OBJETO GetNumberOfPoints]
    ATT CopyComponent 0 [ [OBJETO GetPointData] GetArray $att] $c

    [OBJETO GetPointData] AddArray          ATT
    [OBJETO GetPointData] SetActiveScalars  __actual__
    ATT Delete
  }

  MAPPER SetScalarModeToUsePointData
  MAPPER ScalarVisibilityOn
  eval MAPPER SetScalarRange [[ [OBJETO GetPointData] GetScalars ] GetRange  ]
  VENTANA Render
}

proc CreateCENTERS { } {
  if { [info command CENTERS] == "" } {
    vtkCellCenters  CENTERS
      CENTERS   VertexCellsOff
      CENTERS   SetInput [CROP GetOutput]
  }
  CENTERS   Update
}
proc CreateVECTORS {} {
  if { [info command ONEVECTOR] == "" } {
    vtkArrowSource ONEVECTOR_S
      ONEVECTOR_S SetTipResolution    6
      ONEVECTOR_S SetShaftResolution  1
      ONEVECTOR_S SetShaftRadius      [expr [ONEVECTOR_S GetShaftRadius] * 0.1 ]
    vtkPolyDataNormals ONEVECTOR_N
      ONEVECTOR_N SetInput [ONEVECTOR_S GetOutput]
      ONEVECTOR_N SetFeatureAngle 180
      ONEVECTOR_N Update
    vtkPolyData ONEVECTOR
      ONEVECTOR DeepCopy [ONEVECTOR_N GetOutput]

    ONEVECTOR_S Delete
    ONEVECTOR_N Delete
  }
  if { [info command VECTORS] == "" } {
    vtkGlyph3D  VECTORS
      VECTORS   SetSource ONEVECTOR
      VECTORS   ScalingOn
      VECTORS   OrientOn
      VECTORS   SetScaleModeToScaleByVector
      VECTORS   SetVectorModeToUseVector
      VECTORS   SetColorMode -1
      VECTORS   SetScaleFactor 1
    vtkPolyDataMapper   VECTORSMAPPER
      VECTORSMAPPER    SetInput  [ VECTORS GetOutput ]
    vtkActor            VECTORSACTOR
      VECTORSACTOR     SetMapper VECTORSMAPPER
    RENDERER AddActor VECTORSACTOR
  }
}
proc CreateTENSORS { } {
  if { [info command ONETENSOR] == "" } {
    vtkSphereSource ONETENSOR_S
      ONETENSOR_S SetThetaResolution 8
      ONETENSOR_S SetPhiResolution   8
    vtkPolyDataNormals ONETENSOR_N
      ONETENSOR_N SetInput [ONETENSOR_S GetOutput]
      ONETENSOR_N SetFeatureAngle 180
      ONETENSOR_N Update
    vtkPolyData ONETENSOR
      ONETENSOR DeepCopy [ONETENSOR_N GetOutput]

    ONETENSOR_S Delete
    ONETENSOR_N Delete
  }
  if { [info command TENSORS] == "" } {
    vtkTensorGlyph  TENSORS
      TENSORS SetSource ONETENSOR
      TENSORS ScalingOn
      TENSORS SetScaleFactor 1
      TENSORS SetColorModeToEigenvalues
      TENSORS ColorGlyphsOn
      TENSORS ExtractEigenvaluesOff

    vtkPolyDataMapper TENSORSMAPPER
      TENSORSMAPPER SetInput  [TENSORS GetOutput]
    vtkActor          TENSORSACTOR
      TENSORSACTOR  SetMapper TENSORSMAPPER
    RENDERER AddActor TENSORSACTOR
  }
}
proc DrawTensorsOnNodes { att } {
  MAPPER ScalarVisibilityOff
  if { [info command VECTORSACTOR] != "" } { VECTORSACTOR VisibilityOff  }

  [OBJETO GetPointData] SetActiveTensors $att

  CreateTENSORS
  TENSORS SetInput  [CROP GetOutput]
  TENSORS Update
  TENSORSACTOR VisibilityOn
  VENTANA Render
}

proc DrawTensorsOnCells { att } {
  MAPPER ScalarVisibilityOff
  if { [info command VECTORSACTOR] != "" } { VECTORSACTOR VisibilityOff  }

  CreateCENTERS
  [ [CENTERS GetOutput] GetPointData] SetTensors [[[CROP GetOutput] GetCellData] GetArray $att]

  CreateTENSORS
  TENSORS SetInput  [CENTERS GetOutput]
  TENSORS Update
  TENSORSACTOR VisibilityOn
  VENTANA Render
}

proc DrawVectorsOnNodes { att } {
  MAPPER ScalarVisibilityOff
  if { [info command TENSORSACTOR] != "" } { TENSORSACTOR VisibilityOff  }

  [OBJETO GetPointData] SetActiveVectors $att

  CreateVECTORS
  VECTORS   SetInput [CROP GetOutput]
  VECTORS   Update
  VECTORSACTOR VisibilityOn
  VENTANA Render
}

proc DrawVectorsOnCells { att } {
  MAPPER ScalarVisibilityOff
  if { [info command TENSORSACTOR] != "" } { TENSORSACTOR VisibilityOff  }

  CreateCENTERS
  [ [CENTERS GetOutput] GetPointData] SetVectors [[OBJETO GetCellData] GetArray $att]

  CreateVECTORS
  VECTORS   SetInput [CENTERS GetOutput]
  VECTORS   Update
  VECTORSACTOR VisibilityOn
  VENTANA Render
}

#console show

[vtkDataReader TYPEREADER] SetFileName [lindex $argv 0]
if { [TYPEREADER IsFilePolyData ] } {
  puts polydata
  vtkPolyDataReader READER
  vtkPolyData OBJETO
}
if { [TYPEREADER IsFileUnstructuredGrid ] } {
  vtkUnstructuredGridReader READER
  vtkUnstructuredGrid       OBJETO
}
if { [TYPEREADER IsFileStructuredPoints ] } {
  vtkStructuredPointsReader READER
  vtkStructuredPoints       OBJETO
}
TYPEREADER Delete

READER SetFileName [lindex $argv 0]
READER Update



OBJETO DeepCopy [READER GetOutput]


vtkRenderer RENDERER

vtkRenderWindow VENTANA
  vtkTkRenderWidget .v -width 600 -height 600 -rw VENTANA
  ::vtk::bind_tk_render_widget .v
[VENTANA GetInteractor] SetInteractorStyle [vtkInteractorStyleTrackballCamera STYLE_CAMERA]
VENTANA   AddRenderer RENDERER
pack .v     -expand true  -fill both


vtkPlanes BOUNDS
#vtkClipDataSet CROP
#  CROP SetInput OBJETO
#  CROP SetClipFunction BOUNDS
#  CROP InsideOutOn
vtkExtractGeometry CROP
  CROP SetInput OBJETO
  CROP SetImplicitFunction BOUNDS
  CROP ExtractBoundaryCellsOn

vtkBoxWidget BOX
  BOX SetInteractor [VENTANA GetInteractor]
  BOX SetInput OBJETO
  BOX SetPlaceFactor 1.1
  BOX PlaceWidget
  BOX AddObserver InteractionEvent { BOX GetPlanes BOUNDS }
#      BOX AddObserver InteractionEvent { CROP SetExtent [[ BOX GetProp3D ] GetBounds]  }
  BOX RotationEnabledOn
#  BOX On
  BOX GetPlanes BOUNDS


vtkDataSetMapper MAPPER
  MAPPER SetInput [CROP GetOutput]
  MAPPER ScalarVisibilityOff
vtkActor ACTOR
  ACTOR SetMapper MAPPER
RENDERER  AddActor ACTOR


vtkPolyData     SL_SEEDS
vtkRungeKutta45 SL_INTEGRATOR
vtkStreamLine   SL_STREAMER
  SL_STREAMER SetInput  OBJETO
  SL_STREAMER SetSource SL_SEEDS
  SL_STREAMER SetIntegrator SL_INTEGRATOR
  SL_STREAMER SetMaximumPropagationTime 100
  SL_STREAMER SetIntegrationStepLength .2
  SL_STREAMER SetStepLength .1
#  SL_STREAMER SetNumberOfThreads 1
  SL_STREAMER SetIntegrationDirectionToIntegrateBothDirections
  SL_STREAMER SpeedScalarsOff
  SL_STREAMER VorticityOff
vtkPolyDataMapper SL_MAPPER
  SL_MAPPER SetInput  [SL_STREAMER GetOutput]
vtkActor SL_ACTOR
  SL_ACTOR SetMapper SL_MAPPER
  SL_ACTOR VisibilityOff
RENDERER AddActor SL_ACTOR

vtkLineWidget SL_WIDGET
  SL_WIDGET SetInput OBJETO
  SL_WIDGET SetAlignToYAxis
  SL_WIDGET PlaceWidget
  SL_WIDGET GetPolyData SL_SEEDS
  SL_WIDGET ClampToBoundsOn
  SL_WIDGET SetResolution 50
  SL_WIDGET SetKeyPressActivationValue L
  SL_WIDGET SetInteractor [VENTANA GetInteractor]
  SL_WIDGET AddObserver StartInteractionEvent BeginInteraction
  SL_WIDGET AddObserver InteractionEvent GenerateStreamlines
proc BeginInteraction    {} { SL_ACTOR VisibilityOn }
proc GenerateStreamlines {} { SL_WIDGET GetPolyData SL_SEEDS }


menubutton .att -text "Attributes" -underline 0 -direction below -menu .att.m -relief raised
  menu .att.m -tearoff 1

pack .att


ReadAttributes

package require eW

  set scale 1
  eEntry  .scale  -labeltext "Scale: " -variable scale -step 0.2 \
      -command { if { [info command VECTORS] != "" } { VECTORS   SetScaleFactor $scale }
                 if { [info command TENSORS] != "" } { TENSORS   SetScaleFactor $scale }
                 VENTANA Render
               }

  pack .scale


  eEntry  .opacity  -labeltext "Opacity: " -step 0.1 -from 0 -to 1 -format "%0.02g" \
      -command { [ACTOR GetProperty] SetOpacity [.opacity getValue];
                 VENTANA Render
               }
  .opacity setValue 1

  pack .opacity











