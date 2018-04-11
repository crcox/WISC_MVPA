[m,ix] = max([metadata.ncol]);
xyz = metadata(ix).coords(2).xyz;
radius = 6;

SL0 = generate_searchlights_0(xyz, radius);
SL1 = generate_searchlights_1(xyz, radius);