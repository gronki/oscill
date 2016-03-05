
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

function rainbow, x
    y = 6. * (x - floor(x)) - 3.
    out = bytarr( [ 3, size(y,/dim) ] )
    z = y(*,*,*)
    out[1,*,*,*] = round(exp(-(z)^2 / 1.9) * 255)
    w2 = 1.4
    out[0,*,*,*] = round(exp(-(z-2.)^2 / w2) * 255) + round(exp(-(z+4.)^2 / w2) * 255)
    out[2,*,*,*] = round(exp(-(z+2.)^2 / w2) * 255) + round(exp(-(z-4.)^2 / w2) * 255)
    return, out
end
