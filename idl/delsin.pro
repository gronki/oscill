    @colorr

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    render = 0
    pref = 'norm'

    m = 1.
    hbar = 1.

    a = 1.
	x0 = 0
    s = 0.01
    p = m * 45


    fps = 24
    nrefl = 10
    refltime = 18
    rang_corr = 1.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    ; set_plot,'X'


    !p.charsize = 1.2


    device, decomp=1


    nstat = round(max([300,4*p*a / (hbar*!pi)]))*4
    npoints = 2*nstat

    x = (dindgen(npoints)/(npoints-1.) - 0.5)*a
    dx = x(1) - x(0)

    y = x * dcomplex(0,0)
    j = dcomplex(0,1)
    f = exp(- 0.5 * ( x - x0 ) ^ 2 / s ^ 2) /  sqrt(s*sqrt(!pi) ) $
        * exp( -j * ( x - x0 ) * p / hbar )

	;f = (1 + cos( (4*x/a + 1) * !pi  )  )* ( x lt 0 )
	;f = f / sqrt(total(abs(f)^2.))

    rang = rang_corr * 0.707 * max(abs(f))

    print, 'NAME ', pref

    print, 'states: ', nstat
    print, 'points: ', npoints

    base = dcomplex(0,0) * fltarr(npoints, nstat)
    spec = dcomplex(0,0) * fltarr(nstat)


    for n=1l, nstat do begin
        base(*,n-1) = sin(!pi*n*(x/a+0.5)) * sqrt(2./a)
        spec(n-1) = total(base(*,n-1)*f)*dx
		q = total(abs(spec(0:n-1))^2.)
        print, n, q
		if abs(q-1.) lt 1e-9 then begin
			nstat = n
			break
		endif
    endfor

	base = base(*,0:nstat-1)
	spec = spec(0:nstat-1)

    n = findgen(nstat)+1
    pn = hbar * n * !pi / a
    en = pn ^ 2. / (2.*m)




    if render then begin ;;;;;;;;;;;;;;;;;;;;
        set_plot, 'Z'
        device, z_buf = 0, decomposed = 1, set_resolution = [1600,900], set_pixel_depth = 24
        mk_dir, 'png'
    endif else begin
        set_plot, 'X'

        window, 3, xs=1600, ys=900
        device, decomp = 1
    endelse

    if render eq 0 then fps=fps/3.

    tscal = a * m / (p > (hbar*!pi/a) )
    tmax = nrefl * tscal
    tbins = nrefl * refltime * fps
    print, 'time bins: ', tbins
    t = findgen(tbins)/(tbins-1.) * tmax
    dt = t(1) - t(0)

    ntr = round((tscal*2)/dt)
    trmarg = round(ntr / 9.)
    mmm = rebin(reform(findgen(trmarg),[1,trmarg]),[npoints,trmarg],/samp)
    mmm = exp( -(mmm-trmarg/2.) ^ 2. / ( 2 * (trmarg/4.)^2 ) )
    mmm[*,0] = 1.

    im0 = fltarr(npoints,ntr+trmarg)
    im1 = fltarr(npoints,ntr+trmarg)

    white = rgb(255,255,255)

    t_start = systime(/sec)
    t_msg = t_start

    for i=0l, n_elements(t)-1 do begin
        ar1 = reform(spec*exp(t(i)*en/hbar*j),[1,nstat])
        ar2 = complex(  $
            rebin(real_part(ar1),[npoints,nstat],/samp),    $
            rebin(imaginary(ar1),[npoints,nstat],/samp))
        y = total(ar2 * base,2)
        !p.multi = [0,2,2]
        amp = abs(y)^2
        ph = (atan(y,/phase)/!pi + 1)/2.

        plot, x, amp,yr=[0,1]*rang^2, thick=2, color=white, ystyle=1, xstyle=1
        plot,x,real_part(y),yr=[-1,1]*rang, thick=2, color=rgb(255,255,255), ystyle=1, xstyle=1
        oplot, x, imaginary(y), color=rgb(255, 80, 33), thick=2


        if render ne 0 then begin
            im0 = shift(temporary(im0),0,-1)
            im1 = shift(temporary(im1),0,-1)

            im0(*,ntr:ntr+trmarg-1) = rebin(amp,[npoints,trmarg],/samp)*mmm
            im1(*,ntr:ntr+trmarg-1) = rebin(ph,[npoints,trmarg],/samp)

            plot_image, rgbflip(bytscl(mono2rgb(im0), min=0, max=rang^2)), color=white,   $
                    /nosq, min=0, max=255,  scale=[dx,dt], origin=[x(0),t(i)-ntr*dt]
            plot_image, rgbflip( 1l * rainbow(im1) * mono2rgb(bytscl(im0,min=0, max=rang^2)) / 255. ),$
                    color=white, /nosq, min=0, max=255,  scale=[dx,dt], origin=[x(0),t(i)-ntr*dt]

            write_png, 'png/' + pref + '-' + string(1000000l+i,f='(I0)') + '.png', tvrd(/true)
        endif else begin
            wait, 1/20.
        endelse


        if (systime(/sec) - t_msg) gt 30. then begin
            t_elapsed = systime(/sec) - t_start
             print,n_elements(t)-1-i, '   left: ',  $
                string(1. * t_elapsed / (i+1) * ( n_elements(t)-1 -i) / 60., f='(f0.1)') + ' min ', $
                ' (' + systime(0,t_start + t_elapsed / (i+1) * n_elements(t)) +')'
            t_msg = systime(/sec)
        endif
    endfor

    print, 'FINISHED ', pref

end
