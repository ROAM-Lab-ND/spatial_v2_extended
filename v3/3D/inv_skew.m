function  out = inv_skew( in )

out = 0.5 * [ in(3,2) - in(2,3);
		in(1,3) - in(3,1);
		in(2,1) - in(1,2) ];
