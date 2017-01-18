;**************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************


begin

outputdirectory = "./regridtest/"

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

  ;print(a)


  latitude = a->lat_1
  longitude = a->lon_1
  topo = a->HSURF

  printVarSummary(topo)
  printVarSummary(latitude)
  printVarSummary(longitude)
end if 

topomod = topo({45.933265645:46.510933923},{6.807583514:7.543522346})
topomod!0 = "latitude"      
topomod!1 = "longitude"
dsizes_t = dimsizes(topomod)
dimty= dsizes_t(0)
dimtx=dsizes_t(1)

beta = new((/dimtx,dimty/),double,"No_FillValue")
betax= fspan(1,100,dimtx)
betay= fspan(1,100, dimty)

;printVarSummary(topomod)
;printVarSummary(betax)
beta!0 = "longitude"      ; attribute coordinates
beta!1 = "latitude"
beta&latitude = betay
beta&longitude = betax

beta(:,:)= transpose(topomod(:,:))
;printVarSummary(beta)


alpha = new((/100,100/),double,"No_FillValue")
alphax= ispan(1,100,1)
alphay= ispan(1,100,1)
;printVarSummary(alpha)
;printVarSummary(alphax)
alpha!0 = "latitude"      ; attribute coordinates
alpha!1 = "longitude"
alpha&latitude = alphay
alpha&longitude = alphax
;printVarSummary(alpha)
dsizes_x = dimsizes(alpha)
dimx= dsizes_x(0)
dimy=dsizes_x(1)
alpha(:,:)=10.0
do i=1,10
	alpha({i},{(i):(dimx-i)})=1*i-1
	alpha({dimy-i},{(i):(dimx-i)})=1*i-1
	alpha({(i):(dimy-i)},{i})=1*i-1
	alpha({(i):(dimy-i)},{dimy-i})=1*i-1
end do

b = new((/1000,1000/),double,"No_FillValue")
by= fspan(4,12,1000)
bx= fspan(43,48,1000)
b!0 = "latitude"      ; attribute coordinates
b!1 = "longitude"
b&latitude = by
b&longitude = bx
;printVarSummary(b)
dsizes_x = dimsizes(b)
dimx= dsizes_x(0)
dimy=dsizes_x(1)

printVarSummary(longitude)
printVarSummary(latitude)
printVarSummary(topo)

;c = linint2(longitude,latitude,topo,False,bx,by,0)
c = rcm2rgrid(latitude,longitude,topo,bx,by,0)

  c!0   ="lat"
  c!1   = "lon"
  c&lat = bx
  c&lon = by
printVarSummary(beta)
printVarSummary(b)
printVarSummary(c)


; ===========================
; create plot
; ===========================
 wks   = gsn_open_wks ("png", outputdirectory+"topo")        ; send graphics to PNG file
 wks2   = gsn_open_wks ("png", outputdirectory+"btest")
 wks3   = gsn_open_wks ("png", outputdirectory+"toporegrided_1000")
 wks4   = gsn_open_wks ("png", "transposedtest")
 wks5   = gsn_open_wks ("png", "etest")

 res                 = True                    ; plot mods desired

 res@tiMainString    = "Test"        ; title
 res@cnLevelSpacingF = 0.2                      ; contour interval

 res@cnFillOn        = True                    ; turn on color
 res@cnFillPalette   = "gui_default"           ; choose color map

 res@gsnContourZeroLineThicknessF = 0.2	       ; doubles thickness of zero contour

 plot = gsn_csm_hov(wks, topo, res)
 plot = gsn_csm_hov(wks3, c, res)
  

end