#!/usr/bin/env python3.6

```
Parse output by index sequnence and return index identity
```


# populate index dictonary with known indices
with open(Index, 'r') as ind:
    for line in ind:
        line = line.strip('\n').split('\t')
        Index_dict[line[4]] = line[3]

