
function rgb, r, g, b
    rgb = ((long(r)>0)<255)*1l  $
        + ((long(g)>0)<255)*256l    $
        + ((long(b)>0)<255) * 256l * 256l
    return, rgb
end

function rgb2, rgbarr
    return, rgb(rgbarr(0,*,*,*),rgbarr(1,*,*,*),rgbarr(2,*,*,*))
end

function rgbflip, im
    return, transpose(im,[1,2,0])
end

function mono2rgb, im
    return, rebin(reform(im,[ 1, size(im,/dim) ]), [ 3, size(im,/dim) ], /samp )
end

function sinthr, x
    return, 0.5*sin(((x < 1.) > (-1.))*!pi/2.)+0.5
end
function tanthr, x
    return, atan(x*!pi)  / !pi + 0.5
end

function rainbow, x
    y = 6. * (x - floor(x)) - 3.
    out = fltarr( [ 3, size(y,/dim) ] )
    z = y(*,*,*)

    gp = [-1.9,0.9, 1.8,0.9]
    rp = [0.2,.9, 2.99,1.3]
    bp = [-2.9,1.2, -0.4,1.1]

    out[1,*,*,*] = sinthr((z-gp(0))/gp(1)) - sinthr((z-gp(2))/gp(3))
    out[0,*,*,*] = sinthr((z-rp(0))/rp(1)) - sinthr((z-rp(2))/rp(3)) + sinthr(-(z-rp(2)+6)/rp(3))
    out[2,*,*,*] = sinthr((z-bp(0))/bp(1)) - sinthr((z-bp(2))/bp(3)) + sinthr((z-bp(0)-6)/bp(1))
    return,  bytscl(out,min=0.,max=1.)
end
