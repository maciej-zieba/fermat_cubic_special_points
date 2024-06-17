ring R = (0, u), (x, y, z), dp;

minpoly = 31 + 36u + 27u2 - 4u3 + 9u4 + u6;

poly F = x3 + y3 + z3;

< "data/6_type_points.txt";

list TangentLines;
for (int i = 1; i <= size(Points6type); i++)
{
    TangentLines =
        TangentLines + list(std(poly((Points6type[i][1] ^ 2) * x + (Points6type[i][2] ^ 2) * y + (Points6type[i][3] ^ 2) * z))[1]);
}

list PointsFromTangents;

for (int i = 1; i <= size(TangentLines) - 1; i++)
{
    for (int j = i + 1; j <= size(TangentLines); j++)
    {
        PointsFromTangents = PointsFromTangents + list(std(ideal(TangentLines[i], TangentLines[j])));
    }
}

// Residual intersection points

for (int i = 1; i <= size(Points6type); i++)
{
    ideal ResidualPointIdeal(i) = sat(ideal(F, TangentLines[i]), point_ideal(Points6type[i]))[1];
}

// 27 points, 9 unique points
list ResidualPointsIdeals = ResidualPointIdeal(1..27);
ResidualPointsIdeals = remove_points_duplicates(ResidualPointsIdeals);

list ResidualPoints;

for (int i = 1; i <= size(ResidualPointsIdeals); i++)
{
    ResidualPoints = ResidualPoints + list(point_from_ideal(ResidualPointsIdeals[i]));
}

// Check that residual points are flex points
for (int i = 1; i <= size(ResidualPoints); i++)
{
    list associatedPointOnElliptic = to_elliptic(ResidualPoints[i]);
    if(!is_zero_point(mul_point(3, associatedPointOnElliptic))){
        print("It is not a flex point"); // This should not be printed
    }
}