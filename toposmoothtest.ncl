;**************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************


begin

;******* Adapt this input *******

outputdirectory = "/home/benoit/Documents/NCL/code/regridtest/test/"
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
 head = readAsciiHead("./zinal_astgmt.ascii", 6)     ; there are six header lines, sepcifies the number of lines skipped at the beginning of the file
 ;print(head)
 xdim = str_split(head(0)," ")
 ydim = str_split(head(1)," ")
 xcor = str_split(head(2)," ")
 ycor = str_split(head(3)," ")
 gridsize = str_split(head(4)," ")

 dem = readAsciiTable("./zinal_astgmt.ascii",stringtoint(xdim(1)),"float",6) ; stringtoint(xdim(1)) represent the number of column in the file

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
alphay= ispan(1,dimdemx,1)
alphax= ispan(1,dimdemy,1)
;printVarSummary(alpha)
;printVarSummary(alphax)
alpha!0 = "latitude"      
alpha!1 = "longitude"
alpha&latitude = alphay
alpha&longitude = alphax
;printVarSummary(alpha)
dsizes_x = dimsizes(alpha)
dimx= dsizes_x(0)
dimy=dsizes_x(1)
;alpha(:,:)=0.0
alpha(:,:)=1.0           ; Part that manages the alpha coefficient for the smoothing
imax=400.
do i=1,imax
	alpha({i},{(i):(dimy-i)})=(1/imax)*i-(1/imax)
	alpha({dimx-i},{(i):(dimy-i)})=(1/imax)*i-(1/imax)
	alpha({(i):(dimx-i)},{i})=(1/imax)*i-(1/imax)
	alpha({(i):(dimx-i)},{dimy-i})=(1/imax)*i-(1/imax)
end do

b = new((/dimdemx,dimdemy/),double,"No_FillValue")
bx= fspan(4,12,dimdemx)
by= fspan(43,48,dimdemy)
b!1 = "latitude"      
b!0 = "longitude"
b&latitude = by
b&longitude = bx
;printVarSummary(b)


;printVarSummary(longitude)
;printVarSummary(latitude)
printVarSummary(topo)

cx= fspan(6.0,9.0,dimdemy)
cy= fspan(45.00,48.00,dimdemx)
;c = linint2(longitude,latitude,topo,False,bx,by,0)
c = rcm2rgrid(latitude({45.00:48.00},{6.0:9.0}),longitude({45.00:48.00},{6.0:9.0}),topo({45.00:48.00},{6.0:9.0}),cy,cx,0)

 
  c!0   ="lat"
  c!1   = "lon"
  c&lat = cy
  c&lon = cx

;ccut = transpose(c)
;ccut!1   ="lat"
;ccut!0   = "lon"

printVarSummary(c)
;cmod = linint2(ccut&lat,ccut&lon,ccut,False,alphay,alphax,0)
;printVarSummary(cmod)
printVarSummary(alpha)
printVarSummary(dem)

csmooth=(1.0d-alpha)*c
demsmooth=(alpha)*dem
topoWRFsmooth=csmooth+demsmooth

print("reallocation done")

;****************************

; plotting parameters

;************************************************


 wks   = gsn_open_wks ("png", outputdirectory+"topo")        ; send graphics to PNG file
 wks2   = gsn_open_wks ("png", outputdirectory+"alpha")
 wks3   = gsn_open_wks ("png", outputdirectory+"toporegrided")
 wks4   = gsn_open_wks ("png", outputdirectory+"demtest")
 wks5   = gsn_open_wks ("png", outputdirectory+"topoWRF_smooth")
 wks6   = gsn_open_wks ("png", outputdirectory+"Cosmo_cut")
 wks7   = gsn_open_wks ("png", outputdirectory+"Dem_cut")
 wks8   = gsn_open_wks ("png", outputdirectory+"WRF_cut")

 res                 = True                    ; plot mods desired

 res@tiMainString    = "Test"        ; title
 ;res@cnLevelSpacingF = 0.2                      ; contour interval

 res@cnFillOn        = True                    ; turn on color
 res@cnFillPalette   = "gui_default"           ; choose color map
 ;res@SmoothingOn     = True                    ; Smooth countour boundaries
 res@cnLevelSelectionMode = "ManualLevels"     ; manually specify contour levels
 res@cnMinLevelValF       = 200.             ; min level
 res@cnMaxLevelValF       = 4500.             ; max level
 res@cnLevelSpacingF      =  200.              ; contour interval
 res@gsnContourZeroLineThicknessF = 0.2        ; doubles thickness of zero contour

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
 
 

 ;plot = gsn_csm_hov(wks, topo, res)
 plot = gsn_csm_hov(wks3, c, res)
 plot = gsn_csm_hov(wks2, alpha, resalpha)
 plot = gsn_csm_hov(wks5, topoWRFsmooth, res)
 plot = gsn_csm_hov(wks4, dem, res)

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

;plot = gsn_csm_xy(wks6,c&lon,c(1500,:),resB)
;plot = gsn_csm_xy(wks7,demx,dem(1500,:),resB)
;plot = gsn_csm_xy(wks8,alphax,topoWRFsmooth(1500,:),resB)


printVarSummary(topoWRFsmooth)
printMinMax(topoWRFsmooth,1)
fmtf= dimdemx + "f9.3"
opt        = True
opt@fout = "/home/benoit/Documents/NCL/code/matrix/WRFmatrix_0"
topoWRFsmooth = topoWRFsmooth(::-1,:)
write_matrix(topoWRFsmooth,fmtf,opt)

opt1        = True
opt1@fout = "/home/benoit/Documents/NCL/code/matrix/csmoothmatrix"
write_matrix(csmooth,fmtf,opt1)

opt2        = True
opt2@fout = "/home/benoit/Documents/NCL/code/matrix/demsmoothmatrix"
write_matrix(demsmooth,fmtf,opt2)

; trial 
;asciiwrite("trial.txt",topoWRFsmooth)
end
