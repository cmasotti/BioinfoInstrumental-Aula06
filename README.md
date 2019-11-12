# Bioinformática Instrumental - Módulo II
# Identificação de Variantes Genéticas

## Aula06: Análise de Dados de Variantes Germinativas e Somáticase
Até aqui utilizamos as boas práticas de chamada de variantes preconizadas pelo pacote de ferramentas [GATK](https://software.broadinstitute.org/gatk/) (_Genome Analysis Toolkit_). 

Nesta etapa utilizaremos os dados "limpos", pré-processados (após ApplyBQSR), para iniciar a chamada de variantes propriamente dita. Identificaremos SNPs e INDELs de linhagens tumorais e germinativas para produzir um único conjunto de dados no formato [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (_variant call format_). Por fim, filtraremos variantes potencialmente artefatuais com base em escores de qualidade gerados pelo software de chamada de variantes e pelos demais programas até aqui utlizados.
O workflow completo abordado nesta aula prática:   
![calling](pipeline_aula05.png) 

A seguir, executaremos o passo-a-passo da chamada de variantes e refinamento dos dados:

### PASSO 1: ACESSO AO SERVIDOR REMOTO BIOINFO
Para relembrar como fazer a conexão ao servidor remoto via Putty, reveja os passos [aqui](https://github.com/cmasotti/BioinfoInstrumental-Aula05/blob/master/Acesso_servidor_remoto.pdf).

### PASSO 2: ORGANIZAR ESTRUTURA DE DIRETÓRIOS NO SERVIDOR REMOTO
No prompt da linha de comando no servidor remoto, criar diretórios para pré-processamento.
```bash   
aluno30@ea046e981f34:~$ 
