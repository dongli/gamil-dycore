load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_path")) then
  system("echo '[Error]: Argument file_path is needed!'")
  exit
end if

f = addfile(file_path, "r")

num_lat = dimsizes(f->lat)

; bc_total = vhaeC(f->u(0,1:num_lat-2,:), f->v(0,1:num_lat-2,:))
bc_total = vhaeC(f->u(0,:,:), f->v(0,:,:))

; wks = gsn_open_wks("pdf", systemfunc("basename " + file_path + " .nc") + ".ke_spectrum")
wks = gsn_open_wks("pdf", "ke_spectrum")

ke_spectrum = new(num_lat - 2, typeof(bc_total))
do j = 0, dimsizes(ke_spectrum) - 1
  ke_spectrum(j) = sum(bc_total(0,j,0:j)^2 + bc_total(1,j,0:j)^2 + bc_total(2,j,0:j)^2 + bc_total(3,j,0:j)^2) * 0.25
end do

res = True
res@trXAxisType = "LogAxis"
res@trXMinF = 1
res@trXMaxF = 2d2
res@trYAxisType = "LogAxis"
res@trYMinF = 1.0e-8
res@trYMaxF = 1.0e4
res@xyLineThicknessF = 3.0
res@gsnLeftString = f@desc
res@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")

plot = gsn_csm_xy(wks, ispan(1, dimsizes(ke_spectrum) - 1, 1), ke_spectrum(1:), res)
