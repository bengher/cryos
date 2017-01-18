;**************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************


begin

;******* Adapt this input *******

outputdirectory = "./"
plottype = "png"
Filename = "Smooth_"


inputformat = "nc"
inputdirectory = "./"
year = "2016"
month = "01"
day = "01"
hour = "01"
minute = "00"
seconds = "00"
startdatetime = year+""+month+""+day+""+hour

;********************************

; --- Read Data ----------------------------------------;

if (inputformat .eq. "nc")
  print(inputdirectory+"/laf"+startdatetime+".nc")
  a = addfile(inputdirectory+"/laf"+startdatetime+".nc","r")

  print(a)


  latitude = a->y_1
  longitude = a->x_1
  topo = a->HSURF

  printVarSummary(topo)
  printVarSummary(latitude)
  printVarSummary(longitude)
end if 

; Read dem topography
 head = readAsciiHead("./topo_zinal_2.txt", 6)     ; there are six header lines, sepcifies the number of lines skipped at the beginning of the file
 print(head)
 xdim = str_split(head(0)," ")
 ydim = str_split(head(1)," ")
 xcor = str_split(head(2)," ")
 ycor = str_split(head(3)," ")
 gridsize = str_split(head(4)," ")

 dem = readAsciiTable("./topo_zinal_2.txt",stringtoint(xdim(1)),"float",6) ; stringtoint(xdim(1)) represent the number of column in the file

printVarSummary(dem)

 demx = ispan(stringtoint(xcor(1)),stringtoint(xcor(1))+stringtoint(gridsize(1))*stringtoint(xdim(1))-1,stringtoint(gridsize(1)))
 demy = ispan(stringtoint(ycor(1)),stringtoint(ycor(1))+stringtoint(gridsize(1))*stringtoint(ydim(1))-1,stringtoint(gridsize(1)))
 
 dem = dem(::-1,:)       ; reverse grid
 dem!0 = "latitude"      ; attribute coordinates
 dem!1 = "longitude"
 dem&latitude = demy
 dem&longitude = demx


;*****************************

;Define the boundaries

;*****************************

;south
topo1 = topo({45.933265645},{6.807583514:7.543522346})      ; values corresponding to the meters values, converted wih NAVREF
lont1 = longitude({6.807583514:7.543522346})
dem1 = dem({87050},{551050:608050})                         ;limits of the topo_zinal_2.txt domain +/- 1000 m inside the domain
demx1 = demx

;north
topo2 = topo({46.510933923},{6.807583514:7.543522346})
lont2 = lont1
dem2 = dem({151075},{551050:608050})
demx2 = demx1

;west
topo3 = topo({45.933265645:46.510933923},{6.807583514})
lont3 = latitude({45.933265645:46.5})
dem3 = dem({87050:151075},{551050})
demy3 = demy

;east
topo4 = topo({45.933265645:46.510933923},{7.543522346})
lont4 = lont3
dem4 = dem({87050:151075},{608050})
demy4 = demy3


;****************************

;Define the smoothing/interpolation between the 2 topographies inside the domain

;****************************

;Cosmo

topomod = topo({45.933265645:46.510933923},{6.807583514:7.543522346})

;Topo

demmod = dem({87050:151075},{551050:608050})

;Looping for a coefficient alpha from 0 -> 1, 1 = 1000 m into the domain

topoWRF= demmod

do i = 1,4
  if(i .eq. 1) ;SOUTH
    do alpha = 0,10
      beta= alpha*100
      theta=alpha*0.0009
      print(demmod({87050+beta},{551050:608050}))
      print(topomod({45.9333+theta},{6.807583514:7.543522346}))
      print(topoWRF({87050+beta},{551050:608050}))
      topoWRF({87050+beta},{551050:608050})=(10-alpha)*topomod({45.9333+theta},{6.807583514:7.543522346})+alpha*demmod({87050+beta},{551050:608050}) ;0.0009 -> NAVREF 
    end do
  end if
  if(i .eq. 2) ;NORTH
    do alpha = 0,10
      topoWRF({151075-(alpha*100)},{551050:608050})=(10-alpha)*topomod({46.5109-(alpha*0.0009)},{6.807583514:7.543522346})+alpha*demmod({151075-(alpha*100)},{551050:608050}) 
    end do
  end if
  if(i .eq. 3) ;WEST
    do alpha = 0,10
      topoWRF({87050:151075},{551050+(alpha*100)})=(10-alpha)*topomod({45.933265645:46.510933923},{6.8075+(alpha*0.0013)})+alpha*demmod({87050:151075},{551050+(alpha*100)}) ;0.0013 -> NAVREF 
    end do
  end if
  if(i .eq. 4) ;EAST
    do alpha = 0,10
      topoWRF({87050:151075},{608050-(alpha*100)})=(10-alpha)*topomod({45.933265645:46.510933923},{7.5435-(alpha*0.0013)})+alpha*demmod({87050:151075},{608050-(alpha*100)}) 
    end do
  end if
end do

;****************************

; plotting parameters

;************************************************
 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoW")               ; send graphics to PNG file

; left variable
  resB                        = True 

  resB@pmLegendDisplayMode    = "Always"            ; turn on legend
  resB@pmLegendSide           = "Top"               ; Change location of 
  resB@pmLegendParallelPosF   = .1                  ; move units right
  resB@pmLegendOrthogonalPosF = -0.3                ; move units down
  resB@pmLegendWidthF         = 0.1                ; Change width and
  resB@pmLegendHeightF        = 0.1                ; height of legend.
  resB@lgPerimOn              = False               ; turn off box around
  resB@lgLabelFontHeightF     = .03                 ; label font height
  resB@xyExplicitLegendLabels = (/"Cosmo"/)
  resB@vpWidthF               = 0.6
  resB@vpHeightF              = 0.8
  resB@gsnMaximize            = True
  resB@xyLineThicknesses      = 2.                  ; thicker line
  resB@tiXAxisString = " "  ; axis string
  resB@trYMinF = 400
  resB@trYMaxF = 3600
  resB@xyLineColors      = "blue" 
  resB@xyLineThicknessF = 2
  resB@trXMinF = 6.79
  resB@trXMaxF = 7.56


; right variable
  resT                        = True 

  resT@pmLegendDisplayMode    = "Always"            ; turn on legend
  resT@pmLegendSide           = "Top"               ; Change location of 
  resT@pmLegendParallelPosF   = .1                  ; move units right
  resT@pmLegendOrthogonalPosF = -0.3                ; move units down
  resT@pmLegendWidthF         = 0.1                ; Change width and
  resT@pmLegendHeightF        = 0.1                ; height of legend.
  resT@lgPerimOn              = False               ; turn off box around
  resT@lgLabelFontHeightF     = .03                 ; label font height
  resT@xyExplicitLegendLabels = (/"topo_zinal_2"/)                     
  resT@xyDashPatterns         = 0                   ; dashed line for 2nd
  resT@xyLineThicknesses      = 2                   ; thicker line
  resT@tiXAxisString = " "   ; axis string
  resT@xyLineColors      = "red"
  resT@trXMinF = 550050        				   ; plot1
  resT@trXMaxF = 609050					           ; plot1
  resT@trYMinF = 400
  resT@trYMaxF = 3600
  


; Plot South

;wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoS")               ; send graphics to PNG file
;plot = gsn_csm_x2y2(wks,lont1,demx1,topo1,dem1,resB,resT)
  
; Plot North

;wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoN")               ; send graphics to PNG file
;plot = gsn_csm_x2y2(wks,lont2,demx2,topo2,dem2,resB,resT)
  
; Plot West

 ;resB@trYMinF = 400
 ;resB@trYMaxF = 3600
 ;resB@trXMinF = 45.92
 ;resB@trXMaxF = 46.52

 ;resT@trXMinF = 86050
 ;resT@trXMaxF = 152075
 ;resT@trYMinF = 400
 ;resT@trYMaxF = 3600

;wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoW")               ; send graphics to PNG file
;plot = gsn_csm_x2y2(wks,lont3,demy3,topo3,dem3,resB,resT)
  
; Plot East

;wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoE")               ; send graphics to PNG file
;plot = gsn_csm_x2y2(wks,lont4,demy4,topo4,dem4,resB,resT)
  
; Plot DEM

wks2   = gsn_open_wks (plottype,outputdirectory+Filename+"dem")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Topography demmod"             ; add title
 res@cnLevelSpacingF  = 100
 plot = gsn_contour(wks2,dem({87050:151075},{551050:608050}),res)

wks2   = gsn_open_wks (plottype,outputdirectory+Filename+"topo")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Topography topomod"             ; add title
 res@cnLevelSpacingF  = 100
 plot = gsn_contour(wks2,topo({45.933265645:46.510933923},{6.807583514:7.543522346}),res)

wks2   = gsn_open_wks (plottype,outputdirectory+Filename+"topoWRF")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Topography WRF"             ; add title
 res@cnLevelSpacingF  = 100
 plot = gsn_contour(wks2,topoWRF({87050:151075},{551050:608050}),res) 
;*************************

end
