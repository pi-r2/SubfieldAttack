### Subfield attack on VOX ###

`main.sage` contains a demonstration of the attack using sagemath, and the code used to estimate the gate counts for the attack.
Run it with `sage main.sage`.

In `logs/`, the files `VOX_<name>` contain msolve input corresponding to a VOX instance for parameter set <name>, as described in the paper.
`VOX_<name>.log` contains an msolve output for the previous input. 
These files report the running time but also the full detailed log of the F4 algorithm, as well as the full grevlex Gröbner basis.
