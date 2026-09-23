# sph

Simulacao SPH em C++20 com suporte a OpenMP.

## Estrutura

```text
.
|-- include/   # headers
|-- src/       # codigo-fonte
|-- tests/     # testes de consistencia fisica e numerica
|-- docs/      # documentacao tecnica e evidencias
|-- build/     # objetos e binarios por perfil/formato
`-- output/    # arquivos gerados pela simulacao
```

Os arquivos VTK Legacy sao gravados em `output/vtk/`, o VTKHDF temporal em
`output/vtkhdf/` e os diagnosticos em `output/diagnostics/`.

## Build no Ubuntu/Linux

Dependencias:

- `g++`
- `make`
- runtime OpenMP do GCC
- suporte a `C++20` no compilador
- HDF5 para o formato padrao VTKHDF

Em Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential libgomp1 libhdf5-dev
```

Compilar com o perfil padrao `release` (`-O2`):

```bash
make
```

Perfis de build disponiveis:

```bash
make default    # sem -O explicito
make debug      # -O0 -g
make release    # -O2
make o2         # -O2
make o3         # -O3
make ofast      # -Ofast
```

Para benchmarks, os atalhos mais uteis sao:

```bash
make bench  # roda default, O2, O3 e Ofast e appenda os resultados em docs/benchmark.md
```

Os binarios ficam separados por perfil:

```text
build/debug/sph
build/release/sph
build/default/sph
build/o2/sph
build/o3/sph
build/ofast/sph
```

O link `build/<perfil>/sph` aponta para o binario do ultimo formato compilado
naquele perfil. Os binarios reais ficam em `build/<perfil>/<formato>/sph`.

## Formato de saida

VTKHDF e o formato padrao:

```bash
make PROFILE=o2 OUTPUT_FORMAT=vtkhdf all
```

Alternativas:

```bash
make PROFILE=o2 OUTPUT_FORMAT=legacy all
make PROFILE=o2 OUTPUT_FORMAT=both all
```

Executar:

```bash
make run
```

Executar um perfil especifico:

```bash
make PROFILE=o3 run
make PROFILE=ofast run
```

Os parametros de uma execucao curta podem ser alterados pela CLI:

```bash
./build/o2/sph --steps 1000 --frames 20 --quiet
```

Controlar o numero de threads do OpenMP:

```bash
OMP_NUM_THREADS=8 make run
```

## Testes

Executar a suite de consistencia fisica e numerica:

```bash
make PROFILE=o2 OUTPUT_FORMAT=legacy test
```

Ela verifica geometria das paredes, cobertura do kernel, equilibrio
hidrostatico, continuidade de Adami, pressao por par e Velocity-Verlet.

Limpar binario e objetos:

```bash
make clean
```

Limpar arquivos gerados pela simulacao:

```bash
make clean-output
```

Os diagnosticos principais sao:

```text
output/diagnostics/hydrostatic_initial.csv
output/diagnostics/simulation.csv
```
