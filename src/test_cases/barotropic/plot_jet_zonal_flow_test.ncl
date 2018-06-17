load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  system("echo '[Error]: Argument file_prefix should be set!'")
  exit
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

wks = gsn_open_wks("pdf", file_prefix)

res_common = True
res_common@gsnDraw = False
res_common@gsnFrame = False
res_common@mpCenterLonF = 90
res_common@mpFillOn = False
res_common@mpOutlineOn = False
res_common@pmTickMarkDisplayMode = "Always"
res_common@cnLineLabelsOn = False

res_vor = res_common
res_vor@cnLevelSelectionMode = "ManualLevels"
res_vor@cnMinLevelValF = -0.001
res_vor@cnMaxLevelValF = 0.001
res_vor@cnLevelSpacingF = 2.0e-5
res_vor@gsnContourNegLineDashPattern = 11
res_vor@gsnContourZeroLineThicknessF = 0

res_h = res_common

res_pannel = True
res_pannel@gsnPanelRowSpec = True

f = addfile(file_names(0), "r")
gd0 = f->gd

; Minus perturbation.
; g = 9.80616
; pi = atan(1.0) * 4.0
; ghd = g * 120
; lat2 = pi / 4.0
; alpha = 1.0 / 3.0
; beta = 1.0 / 15.0
; do j = 0, dimsizes(f->lat) - 1
;   do i = 0, dimsizes(f->lon) - 1
;     gd0(0,j,i) = gd0(0,j,i) - ghd * cos(f->lat(j)) * \
;       exp(-((f->lon(i) - pi)  / alpha)^2) * \
;       exp(-((lat2 - f->lat(j)) / beta)^2)
;   end do
; end do 


do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  
  res_vor@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  plot0 = gsn_csm_contour_map(wks, f->vor(0,:,:), res_vor)

  h = (f->gd - gd0) / 9.80616
  copy_VarCoords(f->gd, h)
  h@long_name = "height"
  plot1 = gsn_csm_contour_map(wks, h(0,:,:), res_h)

  gsn_panel(wks, (/plot0,plot1/), (/1,1/), res_pannel)
end do
