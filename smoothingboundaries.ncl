; Smoothing WRF topography

; 2016, Franziska Gerber

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


smoothing = "On" ; else "Off"



;********************************

; --- Read Data ----------------------------------------;

; Read COSMO topography
if (inputformat .eq. "grb1")
 print(inputdirectory+"/laf"+startdatetime+".grb1")
 a = addfile(inputdirectory+"/laf"+startdatetime+".grb1","r")

print(a)

T = a->TMP_GDS0_ISBL
RH = a->R_H_GDS0_ISBL
U = a->U_GRD_GDS0_ISBL
V = a->V_GRD_GDS0_ISBL
gp = a->GP_GDS0_ISBL
latitude = a->g0_lat_1
longitude = a->g0_lon_2
P = a->lv_ISBL0

tk = T(:,{lat},{lon})
tdc = (tk*(1-((tk*log(RH(:,{lat},{lon})/100))/5423))^(-1))-273.15
tc = tk-273.15
z = gp(:,{lat},{lon})
u = U(:,{lat},{lon})
v = V(:,{lat},{lon})
p = P
end if

if (inputformat .eq. "nc")
print(inputdirectory+"/laf"+startdatetime+".nc")
a = addfile(inputdirectory+"/laf"+startdatetime+".nc","r")

print(a)

;T = a->T
;RH = a->RELHUM
;U = a->U
;V = a->V
;gp = a->FI
latitude = a->y_1
longitude = a->x_1
;P = a->P
topo = a->HSURF

printVarSummary(topo)
printVarSummary(latitude)
printVarSummary(longitude)
end if 

; Read dem topography
 head = readAsciiHead("./topo_zinal_2.txt", 6)     ; there are six header lines
 print(head)
 xdim = str_split(head(0)," ")
 ydim = str_split(head(1)," ")
 xcor = str_split(head(2)," ")
 ycor = str_split(head(3)," ")
 gridsize = str_split(head(4)," ")

 dem = readAsciiTable("./topo_zinal_2.txt",stringtoint(xdim(1)),"float",6)
printVarSummary(dem)
 demx = ispan(stringtoint(xcor(1)),stringtoint(xcor(1))+stringtoint(gridsize(1))*stringtoint(xdim(1))-1,stringtoint(gridsize(1)))
 demy = ispan(stringtoint(ycor(1)),stringtoint(ycor(1))+stringtoint(gridsize(1))*stringtoint(ydim(1))-1,stringtoint(gridsize(1)))
 
 dem = dem(::-1,:)       ; reverse grid
 dem!0 = "latitude"      ; attribute coordinates
 dem!1 = "longitude"
 dem&latitude = demy
 dem&longitude = demx


;Read wrf output topography
print(inputdirectory+"/wrfout_d01_2014-07-01_00:00:00")
b = addfile(inputdirectory+"/wrfout_d01_2014-07-01_00:00:00","r")

print(b)

latWRFin = b->XLAT_U
lonWRFin = b->XLONG_V
topoWRF = b->HGT

latWRF = latWRFin(0,:,0)
lonWRF = lonWRFin(0,0,:)

topoWRF&south_north = latWRF
topoWRF&west_east = lonWRF

printVarSummary(topoWRF)
print(latWRF)
print(lonWRF)

;south
;topo1 = topo({45.8},{9.18:10.42})
topo1 = topo({45.8},{6.7:7.6})
lont1 = longitude({6.7:7.6})
;dem1 = dem({144600},{734000:829975})
dem1 = dem({86050},{550050:609050})
demx1 = demx
;topoWRF1 = topoWRF(0,{min(latWRF)},:)
topoWRF1 = topoWRF(0,{min(latWRF)},:)

;north
topo2 = topo({46.5},{6.7:7.6})
lont2 = lont1
dem2 = dem({152075},{550050:609050})
demx2 = demx1
topoWRF2 = topoWRF(0,{max(latWRF)},:)
;west
topo3 = topo({45.8:46.5},{6.7})
lont3 = latitude({45.8:46.5})
dem3 = dem({86050:152075},{550050})
demy3 = demy
topoWRF3 = topoWRF(0,:,{min(lonWRF)})
;east
topo4 = topo({45.8:46.5},{7.6})
lont4 = lont3
dem4 = dem({86050:152075},{609050})
demy4 = demy3
topoWRF4 = topoWRF(0,:,{max(lonWRF)})

; smooth the boundaries

   pi = 4.*atan(1.0)
   rad = (pi/180.)

minlon = 550050-((6.7948-6.7)/(7.5775-6.7948))*(609050-550050)
maxlon = 609050+((7.6-7.5775)/(7.5775-6.7948))*(609050-550050)
lonsCHcor = fspan(minlon,maxlon,46)

topo1CHcor = topo1
printVarSummary(lonsCHcor)
printVarSummary(topo1CHcor&x_1)
topo1CHcor&x_1 = lonsCHcor
dimtopo = dimsizes(dem1)
dem1_before = dem({86050},{550050:609050})

if (smoothing .eq. "On")
  x = 0 ; initialize
  do while (x .lt. dimtopo(0))
    slope = 45. ; initialize
    kk = -1
    do while(slope .gt. 30.)
      kk = kk+1
      dh = abs(topo1CHcor({demx1(x)})-dem(kk,x))
      dl = abs(demx1(kk)-minlon)
      slope = atan(dh/dl)/rad
    end do
    aa = 0 ; initialize
    do while(aa .lt. kk)
      if(topo1CHcor({demx1(x)})-dem(kk,x) .lt. 0)
        dem(aa,x) = topo1CHcor({demx1(x)})+tan(slope*rad)*abs(demx1(aa)-minlon)
      end if
      if(topo1CHcor({demx1(x)})-dem(kk,x) .gt. 0)
        dem(aa,x) = topo1CHcor({demx1(x)})-tan(slope*rad)*abs(demx1(aa)-minlon) 
      end if									
      aa = aa+1
    end do
    x = x+1
  end do
  
end if

;print(slope)
dem1 = dem({86050},{550050:609050})

;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoW")               ; send graphics to PNG file

; left variable
  resB                        = True 
  resB@vpWidthF               = 0.6
  resB@vpHeightF              = 0.8
  resB@gsnMaximize            = True
  resB@xyLineThicknesses      = 2.                  ; thicker line
  resB@tiXAxisString = " "  ; axis string
  resB@trYMinF = 400
  resB@trYMaxF = 3600
  resB@xyLineColors      = "blue" 
  resB@xyLineThicknessF = 2
  resB@trXMinF = 6.7
  resB@trXMaxF = 7.6


; right variable
  resT                        = True                      
  resT@xyDashPatterns         = 0                   ; dashed line for 2nd
  resT@xyLineThicknesses      = 2                   ; thicker line
  resT@tiXAxisString = " "   ; axis string
  resT@xyLineColors      = "red"
  resT@trXMinF = minlon;min(lonWRF);minlon         ; plot1
  resT@trXMaxF = maxlon;max(lonWRF);maxlon         ; plot1
  resT@trYMinF = 400
  resT@trYMaxF = 3600
  


; Plot South
 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoS_WRF")               ; send graphics to PNG file


print(dimsizes(lont1))
print(dimsizes(demx1))
print(dimsizes(lonWRF))
print(dimsizes(topoWRF1))

;  plot = gsn_csm_x2y2(wks,lont1,demx1,topo1,dem1,resB,resT)
  ;plot = gsn_csm_x2y2(wks,lont1,lonWRF,topo1,topoWRF1,resB,resT)
  plot = gsn_csm_x2y2(wks,lonWRF,demx1,topoWRF1,dem1,resB,resT)
;  plot = gsn_csm_xy(wks,lonWRF,topoWRF1,resB)
  

 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoN")               ; send graphics to PNG file
  ;plot = gsn_csm_x2y2(wks,lont2,demx2,topo2,dem2,resB,resT)
  plot = gsn_csm_x2y2(wks,lont2,lonWRF,topo2,topoWRF2,resB,resT)

 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoS_diff")
  plot = gsn_csm_x2y2(wks,lont1,demx1,topo1,dem1_before,resB,resT)

 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoW")               ; send graphics to PNG file
 
  resB@trYMinF = 400
  resB@trYMaxF = 3600
  resB@trXMinF = 45.8
  resB@trXMaxF = 46.5

  resT@trXMinF = min(latWRF);144574.6590
  resT@trXMaxF = max(latWRF);230065.1689
  resT@trYMinF = 400
  resT@trYMaxF = 3600

  ;plot = gsn_csm_x2y2(wks,lont3,demy3,topo3,dem3,resB,resT)
  plot = gsn_csm_x2y2(wks,lont3,latWRF,topo3,topoWRF3,resB,resT)


 wks   = gsn_open_wks (plottype,outputdirectory+Filename+"topoE")               ; send graphics to PNG file
  ;plot = gsn_csm_x2y2(wks,lont4,demy4,topo4,dem4,resB,resT)
  plot = gsn_csm_x2y2(wks,lont4,latWRF,topo4,topoWRF4,resB,resT)



 wks2   = gsn_open_wks (plottype,outputdirectory+Filename+"dem")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Topography"          ; add title
 res@cnLevelSpacingF  = 100
 plot = gsn_contour(wks2,topo({45.9242:46.466},{6.7948:7.5775}),res)


 xwks = gsn_open_wks("x11","gsun02n")

 resources                       = True
 resources@cnMonoFillPattern     = True     ; Turn solid fill back on.
 resources@cnMonoFillColor       = False    ; Use multiple colors.
 resources@cnLineLabelsOn        = False    ; Turn off line labels.
 resources@cnInfoLabelOn         = False    ; Turn off informational
                                              ; label.
 resources@cnLinesOn             = False    ; Turn off contour lines.

   resources@pmLabelBarDisplayMode = "Always" ; Turn on label bar.
   resources@lbPerimOn             = False    ; Turn off perimeter on
                                              ; label bar.
  
   resources@tiMainString    = topo@long_name
   resources@tiMainFont      = 26
    resources@tiXAxisFont     = 26
   resources@tiYAxisFont     = 26
  
  plot = gsn_contour(xwks,dem,False)   

end
