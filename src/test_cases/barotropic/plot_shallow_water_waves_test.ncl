load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.h0.*.nc | head -n1) && echo ${a/\.h0.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if
if (.not. isdefined("NH")) then
  NH = False
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

wks = gsn_open_wks("pdf", file_prefix + ".fft")

; res = True
; res@gsnDraw = False
; res@gsnFrame = False
; res@mpOutlineOn = False
; res@mpFillOn = False
; res@cnLineLabelsOn = False
; res@cnFillOn = True
; res@cnLinesOn = False
; res@cnFillPalette = "BlWhRe"
; res@lbOrientation = "Vertical"
; if (NH) then
;   res@gsnPolar = "NH"
; end if

; res_h = res
; res_h@cnLevelSelectionMode = "ManualLevels"
; res_h@cnMinLevelValF = -1
; res_h@cnMaxLevelValF =  1
; res_h@cnLevelSpacingF = 0.1
; res_h@gsnLeftString = "Height (m)"
; res_h@lbLabelBarOn = True

; res_u = res
; res_u@cnLevelSelectionMode = "ManualLevels"
; res_u@cnMinLevelValF = -1
; res_u@cnMaxLevelValF =  1
; res_u@cnLevelSpacingF = 0.1
; res_u@gsnLeftString = "U (m/s)"
; res_u@lbLabelBarOn = True

; res_v = res
; res_v@cnLevelSelectionMode = "ManualLevels"
; res_v@cnMinLevelValF = -1
; res_v@cnMaxLevelValF =  1
; res_v@cnLevelSpacingF = 0.1
; res_v@gsnLeftString = "V (m/s)"
; res_v@lbLabelBarOn = True

; days = 0
; do t = 0, dimsizes(file_names) - 1, 13
;   system("echo '=> " + file_names(t) + "'")
;   f = addfile(file_names(t), "r")

;   time = f->time

;   h = f->gd(0,:,:) - 5.0e4
;   h = h / max(abs(h))
;   h!0 = "lat"
;   h!1 = "lon"
;   h&lat = f->lat
;   h&lon = f->lon
;   res_h@gsnRightString = "Days " + days
;   ; res_h@gsnRightString = cd_string(time, "Day %j Hour %H")
;   plot0 = gsn_csm_contour_map(wks, h, res_h)

;   u = f->u(0,:,:) / max(abs(f->u(0,:,:)))
;   u!0 = "lat"
;   u!1 = "lon"
;   u&lat = f->lat
;   u&lon = f->lon
;   res_u@gsnRightString = "Days " + days
;   ; res_u@gsnRightString = cd_string(time, "Day %j Hour %H")
;   plot1 = gsn_csm_contour_map(wks, u, res_u)

;   v = f->v(0,:,:) / max(abs(f->v(0,:,:)))
;   v!0 = "lat"
;   v!1 = "lon"
;   v&lat = f->lat
;   v&lon = f->lon
;   res_v@gsnRightString = "Days " + days
;   ; res_v@gsnRightString = cd_string(time, "Day %j Hour %H")
;   plot2 = gsn_csm_contour_map(wks, v, res_v)

;   gsn_panel(wks, (/plot0, plot1, plot2/), (/3, 1/), False)

;   days = days + 13
; end do

; res_hov = True
; res_hov@cnLevelSelectionMode = "ManualLevels"
; res_hov@cnMinLevelValF = -1
; res_hov@cnMaxLevelValF =  1
; res_hov@cnLevelSpacingF = 0.1
; res_hov@cnFillOn = True
; res_hov@cnLinesOn = False
; res_hov@cnLineLabelsOn = False
; ; res_hov@cnFillMode = "RasterFill"
; res_hov@tiXAxisString = "Longitude"
; res_hov@tiYAxisString = "Time (Days)"
; res_hov@gsnLeftString = ""
; res_hov@gsnRightString = ""
; res_hov@vpWidthF = 0.8
; res_hov@vpHeightF = 0.6

; f = addfile(file_names(0), "r")
; start_time_idx = 1142
; end_time_idx = 1202
; v = new((/end_time_idx - start_time_idx, dimsizes(f->lon({0:180}))/), typeof(f->v))
; do t = 0, end_time_idx - start_time_idx - 1
;   system("echo '=> " + file_names(t + start_time_idx) + "'")
;   f = addfile(file_names(t), "r")
;   v(t,:) = f->v(0,{44},{0:180})
;   v(t,:) = (/v(t,:) / max(abs(v(t,:)))/)
; end do
; v!0 = "time"
; v!1 = "lon"
; v&lon = f->lon({0:180})
; plot = gsn_csm_contour(wks, v, res_hov)
; draw(plot)
; frame(wks)

f = addfile(file_names(0), "r")
start_time_idx = 0
end_time_idx = dimsizes(file_names)
v = new((/end_time_idx - start_time_idx/), typeof(f->v))
do t = 0, end_time_idx - start_time_idx - 1
  system("echo '=> " + file_names(t + start_time_idx) + "'")
  f = addfile(file_names(t), "r")
  v(t) = f->v(0,{44},{45})
end do
frq = ezfftf(v)

print(maxind(frq(0,:)))
print((dimsizes(file_names) - 1) / int2flt(maxind(frq(0,:))))
C0 = -0.6061436374682E-06

res_fft = True
res_fft@gsnMaximize = True
res_fft@xyLineThicknessF = 3
res_fft@vpWidthF = 0.8
res_fft@vpHeightF = 0.4

plot = gsn_csm_xy(wks, ispan(0, frq@npts / 2 - 1, 1), frq(0,:), res_fft)