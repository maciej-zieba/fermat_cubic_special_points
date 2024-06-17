// 9-type points on the Fermat cubic

// Division polynomials
ring DP = 0, (x, y, A, B), lp;

poly f(0) = 0;
poly f(1) = 1;
poly f(2) = 2y;
poly f(3) = 3x4 + 6Ax2 + 12Bx - A2;
poly f(4) = 4y * (x6 + 5Ax4 + 20Bx3 - 5A2x2 - 4ABx - 8B2 - A3);

poly f(5) = f(4) * f(2) ^ 3 - f(1) * f(3) ^ 3;
poly f(6) = f(3) * (f(5) * f(2) ^ 2 - f(1) * f(4) ^ 2) / (2y);
poly f(7) = f(5) * f(3) ^ 3 - f(2) * f(4) ^ 3;
poly f(8) = f(4) * (f(6) * f(3) ^ 2 - f(2) * f(5) ^ 2) / (2y);
poly f(9) = f(6) * f(4) ^ 3 - f(3) * f(5) ^ 3;

poly G9 = subst(f(9), A, 0, B, -432);

// Same polynomial, but with substitution y2=x3-432:
poly H9 = 59049x40 
        - 1020366720x37 - 221856x34 * (x3 - 432) ^ 2 
        + 7934371614720x34 + 4762976256x31 * (x3 - 432) ^ 2 
        - 36561584400629760x31 + 261120x28 * (x3 - 432) ^ 4 
        - 40362238476288x28 * (x3 - 432) ^ 2 + 110562231227504394240x28
        - 6465650688x25 * (x3 - 432) ^ 4 + 185984654292025344x25 * (x3 - 432) ^ 2 
        - 229261842673353111896064x25 - 98304x22 * (x3 - 432) ^ 6 
        + 58174803542016x22 * (x3 - 432) ^ 4 - 511267887849538584576x22 * (x3 - 432) ^ 2 
        + 330137053449628481130332160x22 + 2717908992x19 * (x3 - 432) ^ 6 
        - 238190651441676288x19 * (x3 - 432) ^ 4 + 819703589774069034123264x19 * (x3 - 432) ^ 2 
        - 325986759063404580224693698560x19 - 25977774145536x16 * (x3 - 432) ^ 6 
        + 470325011933539860480x16 * (x3 - 432) ^ 4 - 635131350428788308540653568x16 * (x3 - 432) ^ 2 
        + 211239419873086167985601516666880x16 + 93076163257171968x13 * (x3 - 432) ^ 6 
        - 340920657106162509938688x13 * (x3 - 432) ^ 4 + 163094406775023265976766431232x13 * (x3 - 432) ^ 2 
        - 81115937231265088506470982400081920x13 - 64202770792587460608x10 * (x3 - 432) ^ 6 
        - 28339153864171736103124992x10 * (x3 - 432) ^ 4 - 39345215592679823705275952529408x10 * (x3 - 432) ^ 2 
        + 14016833953562607293918185758734155776x10 - 49980871012648564555776x7 * (x3 - 432) ^ 6 
        + 26286957475254572460124667904x7 * (x3 - 432) ^ 4 - 15989932249012271676325123797811200x7 * (x3 - 432) ^ 2 
        - 9487278061310018435678208x4 * (x3 - 432) ^ 6 + 2564250303392575068368290185216x4 * (x3 - 432) ^ 4 
        - 250669387372556915483371507751583744x4 * (x3 - 432) ^ 2 - 565310913446334891615584256x * (x3 - 432) ^ 6 
        - 105500583911008802812866796191744x * (x3 - 432) ^ 4 - 19688940971808106816148452972488032256x * (x3 - 432) ^ 2;

// 9-type points
ring R9 = (0, u), x, dp;
minpoly = u18 - 15u15 + 177u12 - 578u9 + 6747u6 + 642u3 + 343;
poly H9 = imap(DP, H9);
list FH9 = factorize(H9);

// Calculate (x,y) coordinates
list XCoords;
for (int i = 1; i <= deg(H9); i++) {
    XCoords = XCoords + list(-subst(FH9[1][i + 1], x, 0));
}

list YCoords;
for (int i = 1; i <= deg(H9); i++) {
    FH9 = factorize(x^2 - XCoords[i] ^ 3 + 432);
    YCoords = YCoords + list(-subst(FH9[1][2], x, 0), -subst(FH9[1][3], x, 0));
}

// Return to the original form
list Coords;
for (int i = 1; i <= 2 * deg(H9); i++) {
    Coords = Coords + list((36 + YCoords[i]) / 72, (36 - YCoords[i]) / 72, -XCoords[(i + 1) div 2] / 12);
}

list Points9;

for(int i = 1; i <= size(Coords) div 3; i++) {
    Points9 = Points9 + list(list(Coords[3 * i - 2], Coords[3 * i - 1], Coords[3 * i]));
}

// Now we have E[9]/{0} list of points, so we just need to eliminate 3-type points.

list Points9type;

for(int i = 1; i<= size(Points9); i++) {
    list associatedPointOnElliptic = to_elliptic(Points9[i]);
    // Delete all type-3 points
    if(!is_zero_point(mul_point(3, associatedPointOnElliptic))){
        Points9type = Points9type + list(Points9[i]);
    }
}

write(":w data/9_type_points.txt", "list Points9type = ");
for (int i = 1; i <= size(Points9type) - 1; i++) {
    write("data/9_type_points.txt", "list(" + string(Points9type[i][1]) + ", " + string(Points9type[i][2]) + ", " + string(Points9type[i][3]) + "),");
}
write("data/9_type_points.txt", "list(" + string(Points9type[size(Points9type)][1]) + ", " + string(Points9type[size(Points9type)][2]) + ", " + string(Points9type[size(Points9type)][3]) + ");");

ring S = 0, (s, t, u), dp;
list Points9typeST = imap(R9, Points9type);

// List of 72 points in terms of s = e^((2ipi)/9) and t = 3^(1/3)
ideal I = s6 + s3 + 1, t3 - 3;
I = std(I);
list Points9typeSimplified = Points9typeST;

for(int i = 1; i <= size(Points9typeST); i++) {
    Points9typeSimplified[i][1] = reduce(subst(Points9typeST[i][1], u, s + t), I);
    Points9typeSimplified[i][2] = reduce(subst(Points9typeST[i][2], u, s + t), I);
    Points9typeSimplified[i][3] = reduce(subst(Points9typeST[i][3], u, s + t), I);
}

write(":w data/9_type_points_simplified.txt", "list Points9typeSimplified = ");
for (int i = 1; i <= size(Points9typeSimplified) - 1; i++) {
    write("data/9_type_points_simplified.txt", "list(" + string(Points9typeSimplified[i][1]) + ", " + string(Points9typeSimplified[i][2]) + ", " + string(Points9typeSimplified[i][3]) + "),");
}
write("data/9_type_points_simplified.txt", "list(" + string(Points9typeSimplified[size(Points9typeSimplified)][1]) + ", " + string(Points9typeSimplified[size(Points9typeSimplified)][2]) + ", " + string(Points9typeSimplified[size(Points9typeSimplified)][3]) + ");");