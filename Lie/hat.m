function H = hat( h )
  H=zeros(3);
  H(1, 2) = -h(3);
  H(2, 1) =  h(3);
  H(1, 3) =  h(2);
  H(3, 1) = -h(2);
  H(2, 3) = -h(1);
  H(3, 2) =  h(1);
end
