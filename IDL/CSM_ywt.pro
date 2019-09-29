PRO CSM_test
; select the ref and img path
;ref_path = DIALOG_PICKFILE(PATH='/Users/lechen/Desktop')
;img_path = DIALOG_PICKFILE(PATH='/Users/lechen/Desktop')
;ref = read_image(ref_path)
;img = read_image(img_path)
ref = read_image("/Users/lechen/Desktop/WT-1A_PMS_20190618065326_0002695001_001_0001_004_L1/1.tif")
;img = read_image("/Users/lechen/Desktop/MATLAB PROGRAM/match ref/Fig0734(a).tif")

ref_size=size(ref)
;img_size=size(img)
rx = ref_size[1]
ry = ref_size[2]
;ix = img_size[1]
;iy = img_size[2]
Q1 = CSM_gene(rx)
Q2 = CSM_gene(ry)
;img = (transpose(Q1)^(0))#ref#((Q2)^(0))
img = CSM_shift(ref,5,2)
ref =  ref(0:rx-1,0:ry-1)
img =  img(0:rx-1,0:ry-1)
start = systime(/second)
a = ywt_gxr_CSM(ref,img)
stop = systime(/second)
print, stop-start

END

FUNCTION CSM_gene, n
  ; generate the cyclic shift matrix and n is the size of the matrix
  IF n EQ 0 then begin
    print,"the size of n is to small'
    return,0
  ENDIF
  Q = dblarr(1,n)+1
  Q =  diag_matrix(Q)
  Q[0:n-1,0:n-2]=Q[0:n-1,1:n-1]
  Q[0,n-1]=1
  return, Q
END

FUNCTION DFT_ywt, n
  ; generate the DFT matrix and n is the size of the matrix
  res = DCINDGEN(n,n)
  for i = 0, n-1 DO begin
    for j=0, n-1 DO begin
      res[i,j]=dcomplex(cos(2*!Dpi*i*j/n)/sqrt(n),sin(2*!Dpi*i*j/n)/sqrt(n))
    endfor
  endfor
  return, res

END

FUNCTION FS2_ywt, img
  ;my own vision of fftshift in matlab, which is designed for 2-D image
  ;img, the image to be shifted
  img_res = img
  img_size=size(img)
  x_size = img_size[1]
  y_size = img_size[2]
  x_p =0
  y_p =0
  If (x_size mod 2) EQ 1 then begin
    x_p = (x_size+1)/2
  endif ELSE begin
    x_p = x_size/2
  endelse

  If (y_size mod 2) EQ 1 then begin
    y_p = (y_size+1)/2
  endif else begin
    y_p = y_size/2
  endelse

  p1 = img(0:x_p-1,0:y_p-1)
  p2 = img(x_p:x_size-1,0:y_p-1)
  p3 = img(0:x_p-1,y_p:y_size-1)
  p4 = img(x_p:x_size-1,y_p:y_size-1)

  img_res(0:x_size-x_p-1,0:y_size-y_p-1) = p4
  img_res(x_size-x_p:x_size-1,y_size-y_p:y_size-1)=p1
  img_res(0:x_size-x_p-1,y_size-y_p:y_size-1)=p2
  img_res(x_size-x_p:x_size-1,0:y_size-y_p-1)=p3
  return, img_res
END


FUNCTION fft_ywt,img
  img_size=size(img)
  x_size = img_size[1]
  F = DFT_ywt(x_size)
  res = F#img
  return, res
END

FUNCTION ifft_ywt,img
  img_size=size(img)
  x_size = img_size[1]
  F = DFT_ywt(x_size)
  F = conj(transpose(F))
  res = F#img
  return,res
END

FUNCTION ifx_ywt,img
  res = FS2_ywt(fft_ywt(FS2_ywt(img)))
  return, res
END

FUNCTION iifx_ywt,img
  res = FS2_ywt(ifft_ywt(FS2_ywt(img)))
  return, res
END

FUNCTION ify_ywt,img
  temp = transpose(img)
  res = FS2_ywt(fft_ywt(FS2_ywt(temp)))
  res = transpose(res)
  return, res
END

FUNCTION iify_ywt,img
  temp = transpose(img)
  res = FS2_ywt(ifft_ywt(FS2_ywt(temp)))
  res = transpose(res)
  return, res
END



FUNCTION ywt_gxr_CSM, ref, img, matrix=matrix,itenum = itenum;, step=step, subsize=subsize, itenum = itenum
; ref:the reference image 
; img:the image to be matehed
; matrix: the parameter of output pattern
; step: the offset for each calculation
; subsize: the size of the subimage used to calculate the local offset

; read the size of two image and parameter setting
if Keyword_set(matrix) then begin
  matrix = matrix
endif else begin
  matrix = 0
endelse
;if Keyword_set(step) then begin
;  step = step
;endif else begin
;  step = 10
;endelse
;if Keyword_set(subsize) then begin
;  subsize = subsize
;endif else begin
;  subsize = 121
;endelse
if Keyword_set(itenum) then begin
  itenum = itenum
endif else begin
  itenum = 3
endelse

ref_size=size(ref)
img_size=size(img)
rx = ref_size[1]
ry = ref_size[2]
ix = img_size[1]
iy = img_size[2]
IF (rx NE ix) or (ry NE iy) then begin
  print,"the size of the iamge pair does not match, please use the correct image pair"
  return,1 
ENDIF
; resize the image to make the size of the image is (2n+1)*(2n+1)
if (rx mod 2) EQ 0 then begin
  ref = ref(0:rx-2,0:ry-1)
  img = img(0:rx-2,0:ry-1)
  rx = rx-1
endif

if (ry mod 2) EQ 0 then begin
  ref = ref(0:rx-1,0:ry-2)
  img = img(0:rx-1,0:ry-2)
  ry = ry-1
endif
; integer calculation
data1 = ref
data2 = img
;AI = FFT(data1)
;BI = FFT(data2)
AI = ifx_ywt(ify_ywt(data1))
BI = ifx_ywt(ify_ywt(data2))
C = AI*conj(BI)
C = iifx_ywt(iify_ywt(C))
mx = max(abs(C),location)
ind = ARRAY_INDICES(abs(C), location)
resx1 = ind[0]-(rx+1)/2+1
resy1 = ind[1]-(ry+1)/2+2
;print, ind
;print, resx1, resy1
data3 = data1
data4 = data2
rxt = rx
ryt = ry
;print, size(data3)
;print, size(data4)
if resx1 GE 0 then begin
  data3 = data3(resx1:rx-1,0:ry-1)
  data4 = data4(0:rx-1-resx1,0:ry-1)
  rxt = rx-resx1
endif
if resx1 LE 0 then begin
  data3 = data3(0:rx-1+resx1,0:ry-1)
  data4 = data4(0-resx1:rx-1,0:ry-1)
  rxt = rx+resx1
endif

if resy1 GE 0 then begin
  data3 = data3(0:rxt-1,resy1:ry-1)
  data4 = data4(0:rxt-1,0:ry-1-resy1)
  ryt = ry-resy1
endif

if resy1 LE 0 then begin
  data3 = data3(0:rxt-1,0:ry-1+resy1)
  data4 = data4(0:rxt-1,0-resy1:ry-1)
  ryt = ry+resy1
endif
; subpixel calculation
ref_size=size(data3)
rx = ref_size[1]
ry = ref_size[2]
; generate the Q matrix
Q1 = csm_gene(rx)
;print, size(Q1)
Q2 = csm_gene(ry)
; FFT transformation 
F1 = DFT_ywt(rx)
F2 = DFT_ywt(ry)
F1H = conj(transpose(F1))
F2H = conj(transpose(F2))
;print,F1H[1:3,1:3]
D1 = F1H#Q1#F1
D2 = F2H#Q2#F2
FA = F1#data3#F2
FB = F1#data4#F2
C = FA/FB
dA = diag_matrix(D1)
dB = diag_matrix(D2)
aa = (real_part(dA))
bb = IMAGINARY(dA)
angle_c = atan(imaginary(C),real_part(C))
angle_A = dA#(dblarr(1,ry)+1)
angleA_u = atan(imaginary(angle_A),real_part(angle_A))
angle_B = dB#(dblarr(1,rx)+1)
angleB_u = atan(imaginary(angle_B),real_part(angle_B))
angleB_u = transpose(angleB_u)
mask = sqrt(angleA_u*angleA_u+angleB_u*angleB_u) LE 0.5*!Dpi 
;print, size(mask) 
angle_x = angleA_u*mask
angle_y = angleB_u*mask
angle_CC = angle_c*mask
angle_x = reform(angle_x,1,rx*ry)
angle_y = reform(angle_y,1,rx*ry)
angle_CC = reform(angle_CC,rx*ry)
angle_res =[angle_x,angle_y]
res = regress(angle_res,angle_CC,SIGMA=sigma, CONST=const)
;print,res
res_all = dblarr(itenum,2)
res_all(0,0:1)=[resx1-res[0],resy1-res[1]]


for i =1 ,itenum-1 do begin
 ; data4n = ((transpose(Q1))^(res[0]))#data4#((Q2)^(res[1]))
  data4n = CSM_shift(data4,res[0],res[1])
  ;data4n = ((transpose(Q1))^(0.5))#data4#((Q2)^(0.5))
  data3n = data3(1:rx-2,1:ry-2)
  data4n = data4n(1:rx-2,1:ry-2)
  resn = as_subcal(data3n,data4n)
  res = res+resn
  res_all(i,0:1)=[resx1-res[0],resy1-res[1]]
endfor
print,res_all
END

function as_subcal, ref, img
  ref_size=size(ref)
  rx = ref_size[1]
  ry = ref_size[2]
  data3 = ref
  data4 = img
  Q1 = csm_gene(rx)
  Q2 = csm_gene(ry)
  F1 = DFT_ywt(rx)
  F2 = DFT_ywt(ry)
  F1H = conj(transpose(F1))
  F2H = conj(transpose(F2))
  D1 = F1H#Q1#F1
  D2 = F2H#Q2#F2
  FA = F1#data3#F2
  FB = F1#data4#F2
  C = FA/FB
  dA = diag_matrix(D1)
  dB = diag_matrix(D2)
  aa = (real_part(dA))
  bb = IMAGINARY(dA)
  angle_c = atan(imaginary(C),real_part(C))
  angle_A = dA#(dblarr(1,ry)+1)
  angleA_u = atan(imaginary(angle_A),real_part(angle_A))
  angle_B = dB#(dblarr(1,rx)+1)
  angleB_u = atan(imaginary(angle_B),real_part(angle_B))
  angleB_u = transpose(angleB_u)
  mask = sqrt(angleA_u*angleA_u+angleB_u*angleB_u) LE 0.5*!Dpi
  angle_x = angleA_u*mask
  angle_y = angleB_u*mask
  angle_CC = angle_c*mask
  angle_x = reform(angle_x,1,rx*ry)
  angle_y = reform(angle_y,1,rx*ry)
  angle_CC = reform(angle_CC,rx*ry)
  angle_res =[angle_x,angle_y]
  res = regress(angle_res,angle_CC,SIGMA=sigma, CONST=const)
  return,res
END

function CSM_shift, img, sx,sy
  ref_size=size(img)
  rx = ref_size[1]
  ry = ref_size[2]
  Q1 = CSM_gene(rx)
  Q2 = CSM_gene(ry)
  F1 = DFT_ywt(rx)
  F2 = DFT_ywt(ry)
  Q1 = transpose(Q1)
  F1H = conj(transpose(F1))
  F2H = conj(transpose(F2))
  D1 = F1H#Q1#F1
  D2 = F2H#Q2#F2
  d1 = diag_matrix(D1)
  d2 = diag_matrix(D2)
  d1 = d1^(sx)
  d2 = d2^(sy)
  D1 = diag_matrix(d1)
  D2 = diag_matrix(d2)
  res = F1#D1#F1H#img#F2#D2#F2H
  return,res
end
