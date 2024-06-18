proc mul_point(int @m, list @L) {
    list P = @L;
    for (int i = 2; i <= @m; i++) {
        P = add_points(P, @L);
    }
    return (P);
}

proc add_points(list @L1, list @L2) {
    if ((size(@L1) == 1) * (@L1[1] == 0)) {
        return (@L2);
    }

    int tangent = (@L1[1] == @L2[1]) * (@L1[2] == @L2[2]); // boolean check

    list L;

    if (tangent) {
        if (poly(2 * @L1[2]) == 0) {
            return (0);
        }
        poly s = poly(3 * @L1[1] ^ 2) / poly(2 * @L1[2]);
        list L = s ^ 2 - 2 * @L1[1], -@L1[2] - s * (s ^ 2 - 3 * @L1[1]);
    } else {
        if (poly(@L2[1] - @L1[1]) == 0) {
            return (0);
        }
        poly s = poly(@L2[2] - @L1[2]) / poly(@L2[1] - @L1[1]);
        list L = s ^ 2 - @L1[1] - @L2[1], -@L1[2] - s * (s ^ 2 - 2 * @L1[1] - @L2[1]);
    }
    return (L);
}

proc is_zero_point(list @P) {
    for (int i = 1; i <= size(@P); i++) {
        if (@P[i] != 0) {
            return (0);
        }
    }
    return (1);
}

proc to_fermat(list @L) {
    // y^2 = x^3 - 432
    if (@L[2] ^ 2 == @L[1] ^ 3 - 432) {
        list L = poly(@L[2] + 36) / poly(72), 1 - poly(@L[2] + 36) / poly(72), poly(-@L[1]) / poly(12);
        return (L);
    } else {
        return (0);
    }
}

proc to_elliptic(list @L) {
    if (@L[1] + @L[2] == 0) // X + Y == 1
    {
        return (0);
    }
    list uniL = poly(@L[1]) / poly(@L[1] + @L[2]), poly(@L[2]) / poly(@L[1] + @L[2]), poly(@L[3]) / poly(@L[1] + @L[2]);
    // x^3+y^3+z^3=0
    if (uniL[1] ^ 3 + uniL[2] ^ 3 + uniL[3] ^ 3 == 0) {
        list L = -12 * uniL[3],
            36 - 72 * uniL[2];
        return (L);
    } else {
        return (0);
    }
}

proc point_ideal(list @P) {
    matrix M[2][3] = @P[1], @P[2], @P[3],
        x, y, z;
    ideal I = wedge(M, 2);

    return (I)
}

proc set_ideal(list @P) {
    ideal I = 1;
    for (int i = 1; i <= size(@P); i++) {
        I = intersect(I, point_ideal(@P[i]));
    }
    // print(I);

    return (std(I))
}

proc remove_points_duplicates(@L) {
    int i, j, count;
    list seen;

    seen = seen + list(@L[1]);

    for (i = 1; i <= size(@L); i++) {
        count = 0;
        for (j = 1; j <= size(seen); j++) {
            if (std(reduce(std(@L[i]), std(seen[j]))) == 0) {
                count++;
            }
        }
        if (count == 0) {
            seen = seen + list(@L[i]);
        }
    }
    return (seen);
}

proc remove_lines_duplicates(@L) {
    int i, j, count;
    list seen;

    seen = seen + list(@L[1]);

    for (i = 1; i <= size(@L); i++) {
        count = 0;
        for (j = 1; j <= size(seen); j++) {
            if (reduce(@L[i], seen[j]) == 0) {
                count++;
            }
        }
        if (count == 0) {
            seen = seen + list(@L[i]);
        }
    }
    return (seen);
}

proc point_from_ideal(ideal @I) {
    list ones = list(1, 0, 0), list(0, 1, 0), list(0, 0, 1);

    for (int i = 1; i <= 3; i++) {
        list AA = line_coeffs(@I[1]) + line_coeffs(@I[2]) + ones[i];
        matrix b[3][1] = 0, 0, 1;
        matrix A[3][3] = AA[1..size(AA)];
        list L = ludecomp(A);
        list Q = lusolve(L[1], L[2], L[3], b);

        if (size(Q) > 0) {
            if (Q[1] == 1) {
                if ((size(Q[3]) == 1) * (size(Q[2]) == 3)) {
                    return (list(Q[2][1, 1], Q[2][2, 1], Q[2][3, 1]));
                }
            }
        }
    }
}

proc line_coeffs(poly @F) {
    list @L = 0, 0, 0;

    list vars = x, y, z;

    for (int i = 0; i <= 2; i++) {
        if (coef(@F, vars[i mod 3 + 1])[1, 1] == vars[i mod 3 + 1]) {
            @L[i mod 3 + 1] = coef(@F, vars[i mod 3 + 1])[2, 1];
            if (size(coef(@F, vars[i mod 3 + 1])) == 4) {
                @F = coef(@F, vars[i mod 3 + 1])[1, 2] * coef(@F, vars[i mod 3 + 1])[2, 2];
                if (coef(@F, vars[(i + 1) mod 3 + 1])[1, 1] == vars[(i + 1) mod 3 + 1]) {
                    @L[(i + 1) mod 3 + 1] = coef(@F, vars[(i + 1) mod 3 + 1])[2, 1];
                    if (size(coef(@F, vars[(i + 1) mod 3 + 1])) == 4) {
                        @F = coef(@F, vars[(i + 1) mod 3 + 1])[1, 2] * coef(@F, vars[(i + 1) mod 3 + 1])[2, 2];
                        if (coef(@F, vars[(i + 2) mod 3 + 1])[1, 1] == vars[(i + 2) mod 3 + 1]) {
                            @L[(i + 2) mod 3 + 1] = coef(@F, vars[(i + 2) mod 3 + 1])[2, 1];
                            return (@L);
                        }
                    } else {
                        return (@L);
                    }
                }
            } else {
                return (@L);
            }
        }
    }
    return (@L);
}

proc equal_points(ideal @P1, ideal @P2) {
    ideal left_containment = reduce(std(@P1), std(@P2));
    ideal right_containment = reduce(std(@P2), std(@P1));
    return (std(left_containment) == 0) * (std(right_containment) == 0);
}

proc equal_varieties(ideal @P1, ideal @P2) {
    ideal left_containment = reduce(std(@P1), std(@P2));
    ideal right_containment = reduce(std(@P2), std(@P1));
    return (std(left_containment) == 0) * (std(right_containment) == 0);
}

proc member(int @x, list @L) {
    for (int i = 1; i <= size(@L); i++) {
        if (@x == @L[i]) {
            return (1);
        }
    }
    return (0);
}

proc indexes_to_cycles(list @indexes, list @permutations) {
    list cycles;
    list visited;

    for (int i = 1; i <= size(@indexes); i++) {
        if (member(@indexes[i], visited) == 0) {
            list cycle = list(@indexes[i]);
            int current_index = @indexes[i];

            while (1) {
                int next_index = @permutations[current_index];
                if (next_index == @indexes[i]) {
                    break;
                }
                cycle = cycle + list(next_index);
                visited = visited + list(next_index);
                current_index = next_index;
            }
            cycles = cycles + list(cycle);
        }
    }
    return (cycles);
}
