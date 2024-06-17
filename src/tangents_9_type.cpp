ring R = (0, u), (x, y, z), dp;

minpoly = u18 - 15u15 + 177u12 - 578u9 + 6747u6 + 642u3 + 343;

poly F = x3 + y3 + z3;

< "data/9_type_points.txt";

list TangentLines;
for (int i = 1; i <= size(Points9type); i++)
{
    TangentLines =
        TangentLines + list(std(poly((Points9type[i][1] ^ 2) * x + (Points9type[i][2] ^ 2) * y + (Points9type[i][3] ^ 2) * z))[1]);
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

for (int i = 1; i <= size(Points9type); i++)
{
    ideal ResidualPointIdeal(i) = sat(ideal(F, TangentLines[i]), point_ideal(Points9type[i]))[1];
}

// 27 points, 9 unique points
list ResidualPointsIdeals = ResidualPointIdeal(1..72);
ResidualPointsIdeals = remove_points_duplicates(ResidualPointsIdeals);

list ResidualPoints;

for (int i = 1; i <= size(ResidualPointsIdeals); i++)
{
    ResidualPoints = ResidualPoints + list(point_from_ideal(ResidualPointsIdeals[i]));
}

// Check that residual points are type-9 points
for (int i = 1; i <= size(ResidualPoints); i++)
{
    list associatedPointOnElliptic = to_elliptic(ResidualPoints[i]);
    if(!is_zero_point(mul_point(9, associatedPointOnElliptic))){
        print("It is not a type-9 point"); // This should not be printed
    }
}

// Now we will show the set of 72 points is divided into 24 subsets of 3 points
list indexes;
for (int i = 1; i <= size(Points9type); i++)
{
    indexes = indexes + list(i);
}
list permutations;
for (i = 1; i <= size(Points9type); i++)
{
    for (j = 1; j <= size(Points9type); j++)
    {
        if (i != j)
        {
            if (equal_points(point_ideal(Points9type[i]), ResidualPointsIdeals[j]))
            {
                permutations = permutations + list(j);
            }
        }
    }
}

list cycles = indexes_to_cycles(indexes, permutations);

for (i = 1; i <= size(cycles); i++)
{
    // cycles[i][1], cycles[i][2], cycles[i][3]; // Uncomment this line to print the cycles
}

write("data/9_type_points_cycles.txt", "list Cycles = ");
for (int i = 1; i <= size(cycles) - 1; i++)
{
    write("data/9_type_points_cycles.txt", "list(" + string(cycles[i][1]) + ", " + string(cycles[i][2]) + ", " + string(cycles[i][3]) + "),");
}
write("data/9_type_points_cycles.txt", "list(" + string(cycles[size(cycles)][1]) + ", " + string(cycles[size(cycles)][2]) + ", " + string(cycles[size(cycles)][3]) + ");");