ring R = (0, u), (x, y, z, a, b, c, d, e, f), dp;

minpoly = u18 - 15u15 + 177u12 - 578u9 + 6747u6 + 642u3 + 343;

poly F = x3 + y3 + z3;

< "data/9_type_points.txt";

poly T = ax2 + by2 + cz2 + dxy + exz + fyz; // conic

// For each point, make first coordinate = 1

for (int i = 1; i <= size(Points9type); i++)
{
    Points9type[i][2] = number(Points9type[i][2]) / number(Points9type[i][1]);
    Points9type[i][3] = number(Points9type[i][3]) / number(Points9type[i][1]);
    Points9type[i][1] = 1;
}

// Create three conditions for each point

F = subst(F, x, 1);
T = subst(T, x, 1);

for (int i = 1; i <= size(Points9type); i++)
{
    poly @F = subst(F, y, y + Points9type[i][2], z, z + Points9type[i][3]);
    poly @T = subst(T, y, y + Points9type[i][2], z, z + Points9type[i][3]);

    poly T0(i) = subst(@T, y, 0, z, 0);
    @T = @T - T0(i);

    @F = @F / subst(subst(@F, y, 0) / z, z, 0);

    @T = @T - subst(subst(@T, y, 0) / z, z, 0) * @F;

    poly T1(i) = subst(subst(@T, z, 0) / y, y, 0);
    @T = @T - T1(i) * y;

    @T = @T - z * @F * subst(subst(@T, y, 0) / z2, z, 0);

    @T = @T - y * @F * subst(@T / (yz), y, 0, z, 0);

    poly T2(i) = subst(subst(@T, z, 0) / y2, y, 0);
}

// Compute conics

list Conics;

for (int i = 1; i <= size(Points9type) - 1; i++)
{
    for (int j = i + 1; j <= size(Points9type); j++)
    {
        ideal I = T0(i), T1(i), T2(i), T0(j), T1(j), T2(j), ax2 + by2 + cz2 + dxy + exz + fyz, f + 1;
        I = std(I);
        if (size(I) > 1)
        {
            Conics = Conics + list(I[size(I)]);
        }
    }
}

// size(Conics); // There are exactly 324 conics

list ConicsThroughPoints;

for (int i = 1; i <= size(Points9type); i++)
{
    list ConicsThroughPoint;
    for (int j = 1; j <= size(Conics); j++)
    {
        if (subst(Conics[j], x, Points9type[i][1], y, Points9type[i][2], z, Points9type[i][3]) == 0)
        {
            ConicsThroughPoint = ConicsThroughPoint + list(j);
        }
    }
    ConicsThroughPoints = ConicsThroughPoints + list(ConicsThroughPoint); // 9 such conics for each point
}

// The conics passing through a given 9-type point on Fermat cubic intersect in 30 different points besides

list OutsidePointsIdeals;
for(int i = 1; i <= size(Points9type); i++){
    list OutsidePointsIdealsPerPoint;
    for(int j = 1; j <= size(ConicsThroughPoints[i]); j++){
        for(int k = j + 1; k <= size(ConicsThroughPoints[i]); k++){
            ideal I = Conics[ConicsThroughPoints[i][j]], Conics[ConicsThroughPoints[i][k]];
            list OutsidePoint = sat(I, point_ideal(Points9type[i]))[1];
            OutsidePointsIdealsPerPoint = OutsidePointsIdealsPerPoint + list(OutsidePoint[1]);
        }
    }
    OutsidePointsIdealsPerPoint = remove_points_duplicates(OutsidePointsIdealsPerPoint);
    // size(OutsidePointsIdealsPerPoint); // 30 per point
    OutsidePointsIdeals = OutsidePointsIdeals + OutsidePointsIdealsPerPoint;

    list HowManyConicsThroughOutsidePointsPerPoint;
    for(int j = 1; j <= size(OutsidePointsIdealsPerPoint); j++){
        list OutsidePoint = point_from_ideal(OutsidePointsIdealsPerPoint[j]);
        int count = 0;
        for(int k = 1; k <= size(ConicsThroughPoints[i]); k++){
            if(subst(Conics[ConicsThroughPoints[i][k]], x, OutsidePoint[1], y, OutsidePoint[2], z, OutsidePoint[3]) == 0){
                count++;
            }
        }
        HowManyConicsThroughOutsidePointsPerPoint = HowManyConicsThroughOutsidePointsPerPoint + list(count);
    }
    int NumOfPointsWhereTwoConicsIntersect = 0;
    int NumOfPointsWhereThreeConicsIntersect = 0;
    for(int j = 1; j <= size(HowManyConicsThroughOutsidePointsPerPoint); j++){
        if(HowManyConicsThroughOutsidePointsPerPoint[j] == 2){
            NumOfPointsWhereTwoConicsIntersect++;
        }
        if(HowManyConicsThroughOutsidePointsPerPoint[j] == 3){
            NumOfPointsWhereThreeConicsIntersect++;
        }
    }
    // NumOfPointsWhereTwoConicsIntersect; // 27
    // NumOfPointsWhereThreeConicsIntersect; // 3
    
}
OutsidePointsIdeals = remove_points_duplicates(OutsidePointsIdeals);
// size(OutsidePointsIdeals); // 2016 for whole set of points

// Now we will do the same for the set of 540 points

list HowManyConicsThroughOutsidePoints;
for(int i = 1; i <= size(OutsidePointsIdeals); i++){
    list OutsidePoint = point_from_ideal(OutsidePointsIdeals[i]);
    int Count = 0;
    for(int j = 1; j <= size(Conics); j++){
        if(subst(Conics[j], x, OutsidePoint[1], y, OutsidePoint[2], z, OutsidePoint[3]) == 0){
            Count++;
        }
    }
    HowManyConicsThroughOutsidePoints = HowManyConicsThroughOutsidePoints + list(Count);
}

int NumOfPointsWhereTwoConicsIntersect = 0;
int NumOfPointsWhereNineConicsIntersect = 0;

for(int i = 1; i <= size(HowManyConicsThroughOutsidePoints); i++){
    if(HowManyConicsThroughOutsidePoints[i] == 2){
        NumOfPointsWhereTwoConicsIntersect++;
    }
    if(HowManyConicsThroughOutsidePoints[i] == 9){
        NumOfPointsWhereNineConicsIntersect++;
    }
}

// NumOfPointsWhereTwoConicsIntersect; // 1944
// NumOfPointsWhereNineConicsIntersect; // 72

// List of 72 points where 9 conics intersect

list PointsWhereNineConicsIntersect;
ideal PointsWhereNineConicsIntersectIdeal = 1;

for(int i = 1; i <= size(HowManyConicsThroughOutsidePoints); i++){
    if(HowManyConicsThroughOutsidePoints[i] == 9){
        PointsWhereNineConicsIntersect = PointsWhereNineConicsIntersect + list(point_from_ideal(OutsidePointsIdeals[i]));
        PointsWhereNineConicsIntersectIdeal = intersect(PointsWhereNineConicsIntersectIdeal, OutsidePointsIdeals[i]);
    }
}

// PointsWhereNineConicsIntersectIdeal[1]; // xyz

write(":w data/72_points_where_nine_conics_intersect.txt", "list PointsWhereNineConicsIntersect = ");
for (int i = 1; i <= size(PointsWhereNineConicsIntersect) - 1; i++) {
    write("data/72_points_where_nine_conics_intersect.txt", "list(" + string(PointsWhereNineConicsIntersect[i][1]) + ", " + string(PointsWhereNineConicsIntersect[i][2]) + ", " + string(PointsWhereNineConicsIntersect[i][3]) + "),");
}
write("data/72_points_where_nine_conics_intersect.txt", "list(" + string(PointsWhereNineConicsIntersect[size(PointsWhereNineConicsIntersect)][1]) + ", " + string(PointsWhereNineConicsIntersect[size(PointsWhereNineConicsIntersect)][2]) + ", " + string(PointsWhereNineConicsIntersect[size(PointsWhereNineConicsIntersect)][3]) + ");");

// List of 72 points in terms of s = e^((2ipi)/9) and t = 3^(1/3)

ring S = 0, (s, t, u), dp;

list PointsWhereNineConicsIntersectST = imap(R, PointsWhereNineConicsIntersect);

ideal I = s3 - 3, t6 + t3 + 1;
I = std(I);

for(int i = 1; i <= size(PointsWhereNineConicsIntersectST); i++){
    PointsWhereNineConicsIntersectST[i][1] = reduce(subst(PointsWhereNineConicsIntersectST[i][1], u, s + t), I);
    PointsWhereNineConicsIntersectST[i][2] = reduce(subst(PointsWhereNineConicsIntersectST[i][2], u, s + t), I);
    PointsWhereNineConicsIntersectST[i][3] = reduce(subst(PointsWhereNineConicsIntersectST[i][3], u, s + t), I);
}

for(int i = 1; i <= size(PointsWhereNineConicsIntersectST); i++){
    /* Uncomment to print the points in the form of s and t */
    // "Point " + string(i) + ":", 
    // PointsWhereNineConicsIntersectST[i][1] + ", " + PointsWhereNineConicsIntersectST[i][2] + ", " + PointsWhereNineConicsIntersectST[i][3];
    // "";
}