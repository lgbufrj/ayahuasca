
# Projeto Ayahuasca

Este documento descreve a estrutura do repositório, convenções de nomeação e instruções básicas para executar os scripts do projeto. É um guia de referência para colaboradores e para reprodução dos resultados.

## Objetivo

O projeto reúne sequências e dados estruturais de enzimas envolvidas na via biossintética estudada (ayahuasca) e fornece scripts reprodutíveis para buscas por homologia, análises e geração de tabelas/figuras para publicação.

## Convenções principais

- Nomes de pastas e arquivos: inglês, letras minúsculas e underscores em vez de espaços (ex.: `reference_sequence.fasta`, `protein_name`).
- Scripts assumem execução a partir da pasta `scripts/` (caminhos relativos configurados dessa forma).
- Mantenha identificadores externos (UniProt, PubChem) em um único local — `scripts/data.py` — para evitar inconsistências.

Exemplos de nomes de proteínas usados no projeto (não é uma lista exaustiva): `tdc`, `asmt` (use estes nomes como chaves/nomes de pastas).

Exemplos de compostos usados como referência: `tryptamine`, `harmine` (apenas exemplos — não listar tudo aqui).

## Estrutura do repositório

- `genome/` — FASTA(s) do(s) genoma(s) e bancos BLAST (por exemplo: `phased/`, `non_phased/`).
- `transcriptome/` — FASTA(s) do transcriptoma.
- `compounds/` — Dados sobre compostos (estruturas, SDF/MOL, metadados).
- `proteins/` — Dados por proteína; cada proteína tem sua subpasta (`proteins/<protein_name>/`).
- `scripts/` — Scripts de análise e utilitários. Execute-os a partir desta pasta.
- `paper/` — Tabelas, figuras e artefatos preparados para artigos.

Estrutura típica dentro de `proteins/<protein_name>/`:

- `reference/`
	- `sequence.fasta` — sequência de referência (aminoácidos) de um organismo modelo.
	- `structure/` — modelos 3D ou arquivos PDB/mmCIF quando disponíveis.
- `results/` — saídas produzidas (alinhamentos, resultados de BLAST, tabelas).
- `annotations/` — previsões de domínios, anotações funcionais, etc.

## Como executar

1) Crie e ative um ambiente virtual Python (recomendado):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt    # se existir
```

2) Execute scripts a partir da pasta `scripts/`:

```bash
cd /home/pedro/Desktop/projects/ayahuasca/scripts
python nome_do_script.py
```

3) Pipeline de alto nível (ordem típica):

- preparar FASTA(s) em `genome/` e `transcriptome/`
- construir bancos BLAST a partir dos FASTA(s)
- executar buscas de homologia (BLAST/DIAMOND) para as proteínas de interesse
- coletar candidatos em `proteins/<nome>/` e rodar alinhamentos/anotações
- guardar ou prever estruturas em `proteins/<nome>/reference/structure`
- gerar tabelas/figuras em `paper/`

Substitua `nome_do_script.py` pelos nomes reais dos scripts em `scripts/`.

## Dados de referência e IDs

- Centralize IDs (UniProt, PubChem, etc.) em `scripts/data.py` — este arquivo deve ser a fonte canônica dos nomes e mapeamentos.
- Ao adicionar novas proteínas/compostos: crie a pasta `proteins/<name>/` e atualize `scripts/data.py`.

## Contribuição

- Inclua pequenos testes ou exemplos de execução para novos scripts.
- Documente novas dependências externas (versões de BLAST, mafft, IQ-TREE etc.) em `documentation.md` ou em um `tools.md` separado.

## Resolução de problemas e observações

- Scripts supõem diretórios relativos a partir de `scripts/`. Se executar de outra pasta, ajuste caminhos ou altere o working directory.
- Se usar binários externos (BLAST, mafft, IQ-TREE), assegure que estejam instalados e disponíveis no `PATH` ou configure caminhos absolutos.

## Próximos passos sugeridos

- Adicionar `requirements.txt` ou `environment.yml` para facilitar a reprodução do ambiente.
- Criar um `README.md` na raiz com um resumo rápido e link para este `documentation.md`.
- (Opcional) Validar programaticamente `scripts/data.py` e gerar checks automatizados para nomes/IDs.

---

Se quiser, eu posso:

- gerar o `README.md` na raiz automaticamente,
- traduzir este documento para outra variante do português, ou
- abrir `scripts/data.py` e checar os nomes/IDs presentes.


