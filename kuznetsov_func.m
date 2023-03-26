%function v = fun(t,Y)

%  x = Y(1);
%  y = Y(2);
%  v(1) = a.*T.*(1-b.*T) - n.*E.*T;
%  v(2) = s - d.*E + p.*E.*(T./(g+T)) - m.*E.*T;
%  v = v';

%end

function v = fun(t,Y)

  x = Y(1);
  y = Y(2);
  v(1) = 13000 - 0.0412.*E + 0.1245.*E.*(T./((2.019.*10^7)+T)) - (3.422.*(10^-10)).*E.*T;
  v(2) = 0.18.*T.*(1-b.*T) - n.*E.*T;
  v = v';

end
