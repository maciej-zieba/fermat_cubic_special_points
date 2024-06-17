ring R9 = (0, u), (x, y, z), dp;
minpoly = u18 - 15u15 + 177u12 - 578u9 + 6747u6 + 642u3 + 343;

< "data/9_type_points.txt";

ideal I = slimgb(set_ideal(Points9type));

poly Q = I[2]; // deq 24

// equal_varieties(ideal(x3 + y3 + z3, Q), I); // true, x3 + y3 + z3 and 24 deg polynomial

// Write 24 deg poly to the file

write(":w data/24_deg_curve.txt", "poly curve24deg = ");
write("data/24_deg_curve.txt", string(I[2]) + ";");