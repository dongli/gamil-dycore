load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  system("echo '[Error]: Argument file_prefix should be set!'")
  exit
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")


wks = gsn_open_wks("pdf", file_prefix)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@mpOutlineOn = False
res@mpFillOn = False
res@cnConstFEnableFill = True
res@gsnLeftString = "Height (m)"

res_panel = True

f0 = addfile(file_names(0), "r")
gd0 = f0->gd(0,:,:)
gd0 = gd0 / 9.80616
do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  gd = f->gd(0,:,:)
  gd = gd / 9.80616
  res@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, gd, res)
  error = gd - gd0
  copy_VarCoords(gd, error)
  plot1 = gsn_csm_contour_map(wks, error, res)
  gsn_panel(wks, (/plot0,plot1/), (/2,1/), res_panel)
end do