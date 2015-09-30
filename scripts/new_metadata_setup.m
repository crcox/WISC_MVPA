metadata0 = metadata;
metadata = metadata0;
for i=1:10
  metadata(i).filter(1).label = 'rowfilter';
  metadata(i).filter(1).dimension = 1;
  metadata(i).filter(1).filter = metadata(i).TRfilter;
  metadata(i).filter(1).note = 'outlying rows';
  metadata(i).filter(2).label = 'colfilter';
  metadata(i).filter(2).dimension = 2;
  metadata(i).filter(2).filter = metadata(i).voxfilter;
  metadata(i).filter(2).note = 'outlying columns';
  metadata(i).filter(3).label = 'faces';
  metadata(i).filter(3).dimension = 1;
  metadata(i).filter(3).filter = metadata(i).TrueFaces;
  metadata(i).filter(3).note = 'Trials with face stimuli';
  metadata(i).filter(4).label = 'places';
  metadata(i).filter(4).dimension = 1;
  metadata(i).filter(4).filter = metadata(i).TruePlaces;
  metadata(i).filter(4).note = 'Trials with place stimuli';
  metadata(i).filter(5).label = 'objects';
  metadata(i).filter(5).dimension = 1;
  metadata(i).filter(5).filter = metadata(i).TrueObjects;
  metadata(i).filter(5).note = 'Trials with object stimuli';
end

for i=1:10
  metadata(i).coords(1).orientation = 'tlrc';
  metadata(i).coords(1).ind = [];
  metadata(i).coords(1).ijk = [];
  metadata(i).coords(1).xyz = metadata(i).xyz_tlrc;
  metadata(i).coords(1).note = 'manually warped tlrc coords.';
end

for i=1:10
  metadata(i).nrow = numel(metadata(i).TrueFaces);
  metadata(i).ncol = size(metadata(i).coords(1).xyz, 1);
end

metadata = rmfield(metadata, 'CVBLOCKS');
metadata = rmfield(metadata, 'xyz_tlrc');
metadata = rmfield(metadata, 'TRfilter');
metadata = rmfield(metadata, 'voxfilter');