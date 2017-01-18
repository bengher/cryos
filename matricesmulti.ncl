;**************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************


begin


alpha = new((/100,100/),double,"No_FillValue")
alphax= ispan(1,100,1)
alphay= ispan(1,100,1)
printVarSummary(alpha)
printVarSummary(alphax)
alpha!0 = "latitude"      ; attribute coordinates
alpha!1 = "longitude"
alpha&latitude = alphay
alpha&longitude = alphax
;printVarSummary(alpha)
dsizes_x = dimsizes(alpha)
dimx= dsizes_x(0)
dimy=dsizes_x(1)
alpha(:,:)=1.0
do i=1,10
	alpha({i},{(i):(dimx-i)})=0.1*i-0.1
	alpha({dimy-i},{(i):(dimx-i)})=0.1*i-0.1
	alpha({(i):(dimy-i)},{i})=0.1*i-0.1
	alpha({(i):(dimy-i)},{dimy-i})=0.1*i-0.1
end do

b = new((/100,100/),double,"No_FillValue")
bx= ispan(1,100,1)
by= ispan(1,100,1)
b!0 = "latitude"      ; attribute coordinates
b!1 = "longitude"
b&latitude = by
b&longitude = bx
;printVarSummary(b)
dsizes_x = dimsizes(b)
dimx= dsizes_x(0)
dimy=dsizes_x(1)
b(:,:)=10.0

c=alpha*b
d=(1-alpha)*b
e=d+c


; ===========================
; create plot
; ===========================
 wks   = gsn_open_wks ("png", "alphatest")        ; send graphics to PNG file
 wks2   = gsn_open_wks ("png", "btest")
 wks3   = gsn_open_wks ("png", "ctest")
 wks4   = gsn_open_wks ("png", "dtest")
 wks5   = gsn_open_wks ("png", "etest")

 res                 = True                    ; plot mods desired

 res@tiMainString    = "Test"        ; title
 res@cnLevelSpacingF = 1.                      ; contour interval

 res@cnFillOn        = True                    ; turn on color
 res@cnFillPalette   = "gui_default"           ; choose color map

 res@gsnContourZeroLineThicknessF = 2.	       ; doubles thickness of zero contour

 plot = gsn_csm_hov(wks, alpha, res)
 plot = gsn_csm_hov(wks2, b, res)
 plot = gsn_csm_hov(wks3, c, res)
 plot = gsn_csm_hov(wks4, d, res)
 plot = gsn_csm_hov(wks5, e, res)


end
