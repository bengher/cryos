;**************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************


begin

;******* Adapt this input *******

outputdirectory = "/home/benoit/Documents/NCL/code/regridtest/"
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

;********************************


if (inputformat .eq. "nc")
  print(inputdirectory+"/laf"+startdatetime+".nc")
  a = addfile(inputdirectory+"/laf"+startdatetime+".nc","r")

  ;print(a)


  latitude = a->lat_1
  longitude = a->lon_1
  topo = a->HSURF

  ;printVarSummary(topo)
  ;printVarSummary(latitude)
  ;printVarSummary(longitude)
end if 

; Read dem topography
 head = readAsciiHead("./topo_zinal_2.txt", 6)     ; there are six header lines, sepcifies the number of lines skipped at the beginning of the file
 ;print(head)
 xdim = str_split(head(0)," ")
 ydim = str_split(head(1)," ")
 xcor = str_split(head(2)," ")
 ycor = str_split(head(3)," ")
 gridsize = str_split(head(4)," ")

 dem = readAsciiTable("./topo_zinal_2.txt",stringtoint(xdim(1)),"float",6) ; stringtoint(xdim(1)) represent the number of column in the file

;printVarSummary(dem)

 demx = ispan(stringtoint(xcor(1)),stringtoint(xcor(1))+stringtoint(gridsize(1))*stringtoint(xdim(1))-1,stringtoint(gridsize(1)))
 demy = ispan(stringtoint(ycor(1)),stringtoint(ycor(1))+stringtoint(gridsize(1))*stringtoint(ydim(1))-1,stringtoint(gridsize(1)))
 
 dem = dem(::-1,:)       ; reverse grid
 dem!0 = "latitude"      ; attribute coordinates
 dem!1 = "longitude"
 dem&latitude = demy
 dem&longitude = demx
 ;printVarSummary(demx)
 dsizes_dem = dimsizes(dem)
 dimdemx= dsizes_dem(0)
 dimdemy=dsizes_dem(1)



;*****************************

;Define the smoothing/interpolation between the 2 topographies inside the domain

;****************************

alpha = new((/dimdemx,dimdemy/),double,"No_FillValue")
alphax= ispan(1,dimdemx,1)
alphay= ispan(1,dimdemy,1)
;printVarSummary(alpha)
;printVarSummary(alphax)
alpha!1 = "latitude"      ; attribute coordinates
alpha!0 = "longitude"
alpha&latitude = alphay
alpha&longitude = alphax
;printVarSummary(alpha)
dsizes_x = dimsizes(alpha)
dimx= dsizes_x(0)
dimy=dsizes_x(1)
alpha(:,:)=1.0
do i=1,100
	alpha({i},{(i):(dimy-i)})=0.01*i-0.01
	alpha({dimx-i},{(i):(dimy-i)})=0.01*i-0.01
	alpha({(i):(dimx-i)},{i})=0.01*i-0.01
	alpha({(i):(dimx-i)},{dimy-i})=0.01*i-0.01
end do

;****************************

; plotting parameters

;************************************************


 wks   = gsn_open_wks ("png", outputdirectory+"topo")        ; send graphics to PNG file
 wks2   = gsn_open_wks ("png", outputdirectory+"alphatest")
 wks3   = gsn_open_wks ("png", outputdirectory+"toporegrided")
 wks4   = gsn_open_wks ("png", outputdirectory+"demtest")
 wks5   = gsn_open_wks ("png", outputdirectory+"alphacut")

resB                        = True 

  resB@pmLegendDisplayMode    = "Always"            ; turn on legend
  resB@pmLegendSide           = "Top"               ; Change location of 
  resB@pmLegendParallelPosF   = .1                  ; move units right
  resB@pmLegendOrthogonalPosF = -0.3                ; move units down
  resB@pmLegendWidthF         = 0.1                ; Change width and
  resB@pmLegendHeightF        = 0.1                ; height of legend.
  resB@lgPerimOn              = False               ; turn off box around
  resB@lgLabelFontHeightF     = .03                 ; label font height
    resB@vpWidthF               = 0.6
  resB@vpHeightF              = 0.8
  resB@gsnMaximize            = True
  resB@xyLineThicknesses      = 2. 


 res                 = True                    ; plot mods desired

 res@tiMainString    = "Test"        ; title
 res@cnLevelSpacingF = 0.2                      ; contour interval

 res@cnFillOn        = True                    ; turn on color
 res@cnFillPalette   = "gui_default"           ; choose color map

 res@gsnContourZeroLineThicknessF = 0.2	       ; doubles thickness of zero contour

 resalpha                 = True                    ; plot mods desired

 resalpha@tiMainString    = "Test"        ; title
 resalpha@cnLevelSpacingF = 0.2                      ; contour interval

 resalpha@cnFillOn        = True                    ; turn on color
 resalpha@cnFillPalette   = "gui_default"           ; choose color map
 resalpha@cnLevelSelectionMode = "ManualLevels"     ; manually specify contour levels
 resalpha@cnMinLevelValF       = 0             ; min level
 resalpha@cnMaxLevelValF       = 1.             ; max level
 resalpha@cnLevelSpacingF      =  0.05

 resalpha@gsnContourZeroLineThicknessF = 0.2        ; doubles thickness of zero contour

 plot = gsn_csm_hov(wks2, alpha, resalpha)
 ;plot = gsn_csm_hov(wks4, dem, res)
 ;plot = gsn_csm_xy(wks5,alphax,alpha(:,1500),resB)


end
