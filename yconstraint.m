function yout = yconstraint(yin, fimg)
  yout = ifft2( abs(fimg) .* exp( 1j * angle(fft2(yin)) ), 'symmetric' );
end
