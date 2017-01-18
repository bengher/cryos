;**************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************


begin

  N    = 5    ; number of rows
  M    = 7    ; number of columns
  ave  = 0.0
  std  = 5.0
  x    = random_normal (ave, std, (/N,M/)) 
  

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


  alpha = new((/dimdemx,dimdemy/),double,"No_FillValue")
alphay= ispan(1,dimdemx,1)
alphax= ispan(1,dimdemy,1)
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
alpha(:,:)=1.0
do i=1,100
	alpha({i},{(i):(dimy-i)})=0.01*i-0.01
	alpha({dimx-i},{(i):(dimy-i)})=0.01*i-0.01
	alpha({(i):(dimx-i)},{i})=0.01*i-0.01
	alpha({(i):(dimx-i)},{dimy-i})=0.01*i-0.01
end do

  printVarSummary(x)
  printVarSummary(alpha)

  fmtf   = dimdemy + "f7.2"
  opt        = True
  opt@fout = "/home/benoit/Documents/NCL/code/matrix/WRFmatrix"
  write_matrix (alpha, fmtf, opt)




end