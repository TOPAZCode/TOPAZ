M4 = Union[ Permutations[{5, 4, 3, 2}], Permutations[{5, 4, 3, 1}], Permutations[{5, 4, 2, 1}], Permutations[{5, 3, 2, 1}], Permutations[{4, 3, 2, 1}] ]

M3a = Union[Permutations[{5, 4, 3}], Permutations[{5, 4, 2}], Permutations[{5, 3, 2}], Permutations[{4, 3, 2}]];
M3b = Union[Permutations[{5, 4, 3}], Permutations[{5, 4, 1}], Permutations[{5, 3, 1}], Permutations[{4, 3, 1}]];
M3c = Union[Permutations[{5, 4, 2}], Permutations[{5, 4, 1}], Permutations[{5, 2, 1}], Permutations[{4, 2, 1}]];
M3d = Union[Permutations[{5, 3, 2}], Permutations[{5, 3, 1}], Permutations[{5, 2, 1}], Permutations[{3, 2, 1}]];
M3e = Union[Permutations[{4, 3, 2}], Permutations[{4, 3, 1}], Permutations[{4, 2, 1}], Permutations[{3, 2, 1}]];
M3 = Union[M3a, M3b, M3c, M3d, M3e]


keyList = {};
OffSet = 37;
Base = 5;
For[i = 1, i <= Length[M4], i++,
  key = M4[[i, 1]]*Base^3 + M4[[i, 2]]*Base^2 + M4[[i, 3]]*Base^1 + M4[[i, 4]]*Base^0 - OffSet;
  keyList = Append[keyList, key];
];
For[i = 1, i <= Length[M3], i++,
  key =  M3[[i, 1]]*Base^2 + M3[[i, 2]]*Base^1 + M3[[i, 3]]*Base^0 - OffSet;
  keyList = Sort[Append[keyList, key]];
];




(* write fortran code *)
$FortranOutput = OpenWrite["/home/schulze/projects/ttbjets/Caching.f90"];
j=1;
WriteString[$FortranOutput,"integer, parameter :: Cache_PartKey(1:"<>ToString[Max[keyList]]<>") = (\\ "];
For[i = 1, i <= Max[keyList], i++,
  If[i == keyList[[j]],
        If[ i != Max[keyList],
            WriteString[$FortranOutput,ToString[j]<>","];
            j++;
        ,
            WriteString[$FortranOutput,ToString[j]<>" \\)"];
         ];
    ,
        WriteString[$FortranOutput,"0,"];
    ];
];

Close[$FortranOutput];
