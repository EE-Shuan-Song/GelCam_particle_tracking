function agg = check_agg_area(area1,area2,area3)

area = [area1 area2 area1+area2 area3];
[B,I] = sort(area);
switch I(4)
    case 4
        agg = 3;
    case 1
        agg = find(I == 2);
    case 2
        agg = 0;
    case 3
        if (B(4) - area3) < (area3 - B(2))
            agg = 3;
        else 
            agg = 0;
        end
end

end