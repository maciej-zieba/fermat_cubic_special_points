ring R = (0, u), (x, y, z, a, b, c, d, e, f), dp;

minpoly = 31 + 36u + 27u2 - 4u3 + 9u4 + u6;

poly F = x3 + y3 + z3;

< "data/6_type_points.txt";

poly T = ax2 + by2 + cz2 + dxy + exz + fyz; // conic

// For each point, make first coordinate = 1

for (int i = 1; i <= size(Points6type); i++)
{
    Points6type[i][2] = number(Points6type[i][2]) / number(Points6type[i][1]);
    Points6type[i][3] = number(Points6type[i][3]) / number(Points6type[i][1]);
    Points6type[i][1] = 1;
}

// Create three conditions for each point

F = subst(F, x, 1);
T = subst(T, x, 1);

for (int i = 1; i <= size(Points6type); i++)
{
    poly @F = subst(F, y, y + Points6type[i][2], z, z + Points6type[i][3]);
    poly @T = subst(T, y, y + Points6type[i][2], z, z + Points6type[i][3]);

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

for (int i = 1; i <= size(Points6type) - 1; i++)
{
    for (int j = i + 1; j <= size(Points6type); j++)
    {
        ideal I = T0(i), T1(i), T2(i), T0(j), T1(j), T2(j), ax2 + by2 + cz2 + dxy + exz + fyz, f + 1;
        I = std(I);
        if (size(I) > 1)
        {
            Conics = Conics + list(I[size(I)]);
        }
    }
}

// size(Conics); // There are exactly 108 conics

list ConicsThroughPoints;

for (int i = 1; i <= size(Points6type); i++)
{
    list ConicsThroughPoint;
    for (int j = 1; j <= size(Conics); j++)
    {
        if (subst(Conics[j], x, Points6type[i][1], y, Points6type[i][2], z, Points6type[i][3]) == 0)
        {
            ConicsThroughPoint = ConicsThroughPoint + list(j);
        }
    }
    ConicsThroughPoints = ConicsThroughPoints + list(ConicsThroughPoint); // 8 such conics for each point
}

// The conics passing through a given sextactic point on Fermat cubic intersect in 24 different points besides

list OutsidePointsIdeals;
for(int i = 1; i <= size(Points6type); i++){
    list OutsidePointsIdealsPerPoint;
    for(int j = 1; j <= size(ConicsThroughPoints[i]); j++){
        for(int k = j + 1; k <= size(ConicsThroughPoints[i]); k++){
            ideal I = Conics[ConicsThroughPoints[i][j]], Conics[ConicsThroughPoints[i][k]];
            list OutsidePoint = sat(I, point_ideal(Points6type[i]))[1];
            OutsidePointsIdealsPerPoint = OutsidePointsIdealsPerPoint + list(OutsidePoint[1]);
        }
    }
    OutsidePointsIdealsPerPoint = remove_points_duplicates(OutsidePointsIdealsPerPoint);
    // size(OutsidePointsIdealsPerPoint); // 24 per point
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
    // NumOfPointsWhereTwoConicsIntersect; // 22
    // NumOfPointsWhereThreeConicsIntersect; // 2
    
}
OutsidePointsIdeals = remove_points_duplicates(OutsidePointsIdeals);
// size(OutsidePointsIdeals); // 540 for whole set of points

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
int NumOfPointsWhereSixConicsIntersect = 0;
int NumOfPointsWhereNineConicsIntersect = 0;

for(int i = 1; i <= size(HowManyConicsThroughOutsidePoints); i++){
    if(HowManyConicsThroughOutsidePoints[i] == 2){
        NumOfPointsWhereTwoConicsIntersect++;
    }
    if(HowManyConicsThroughOutsidePoints[i] == 6){
        NumOfPointsWhereSixConicsIntersect++;
    }
    if(HowManyConicsThroughOutsidePoints[i] == 9){
        NumOfPointsWhereNineConicsIntersect++;
    }
}

// NumOfPointsWhereTwoConicsIntersect; // 486
// NumOfPointsWhereSixConicsIntersect; // 36
// NumOfPointsWhereNineConicsIntersect; // 18

