
@colorr


nx = 2000
a = 20.
hbar = 1.
mass = 1.
omega = 1.


x = (dindgen(nx)/(nx-1.) - 0.5) * a
volt = ((x lt -a/2.) or (x gt a/2.)) * 1d8
volt = mass * omega^2 * x^2
volt =  (x gt 0) * 0.5
volt = (sin( 3 * !pi * (x/a + 1/2.)) * ((x gt -a/2.) and (x lt a/2.)) + 1.)
volt = ((x lt -a/3.) or (x gt a/3.)) * 6.
volt = (x/a + 0.5)*8
volt = (cos( 3*2*!pi * (x/a+1/2.))+1)*3
volt =  16*exp(-x^2 / (0.02*a)^2)

name = 'szpila'

volt(0:2) = 1e6 & volt(nx-3:nx-1) = volt(0:2)
dx = x(2) - x(1)


ham = dblarr(nx,nx)


for j = 0,nx-1 do begin
  y = x * 0.
  y[j] = 1/sqrt(dx)

  p2 = -convol(y,[1,-2,1])/(dx*hbar)^2
  ham[*,j] = (p2 / (2*mass) + volt*y) * sqrt(dx)
endfor

ix = where(ham ne transpose(ham),cw)
if cw ne 0 then ham(ix) = 0
; print,ix
; print,ham(ix)
; print,(transpose(ham))(ix)

e = EIGENQL( ham, EIGENVEC=ev )
ev = ev / sqrt(dx)

set_plot, 'Z'
!P.charsize = 1.4
device, z_buf = 0, decomposed = 1, set_resolution = [1920,1080], set_pixel_depth = 24
mk_dir, 'png2'

!p.multi = [0,4,2]
device,decomp=0
loadct,7

plot, x(4:nx-5), volt(4:nx-5), title='Potencjal'

nlist = [1,2,3,4,6,10,15]
for i = 0,6  do begin
  plot, x, ev[*,nx-nlist(i)], title="Stan "+string(nlist(i),f='(I0)')
  ; oplot, x, sqrt(2./a) * sin( !pi * i * (x/a + 0.5) ), color=200
endfor
for i = 0, nx-1 do print, i, total((ev[*,nx-i-1])^2.)*dx, e(nx-i-1), ((hbar * !pi * (i+1))^2) / (2*mass*a^2)

write_png, 'png2/'+name+'.png',tvrd(/true)
; e = HQR(ELMHES(ham), /DOUBLE)
; vec = eigenvec( ham, e )

end
