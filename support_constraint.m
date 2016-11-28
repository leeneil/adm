function xout = support_constraint(xin, mask)
  xin(~mask) = 0;
  xin(xin < 0) = 0;
  xout = xin;
end
