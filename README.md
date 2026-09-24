# Chaos

Projeto Python.

## Requisitos

- Python 3.10+

## Instalação

```bash
# Crie e ative um ambiente virtual
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate   # Windows

# Instale as dependências
pip install -e ".[dev]"
```

## Uso

```bash
python -m chaos
```

## Desenvolvimento

Rodar os testes:

```bash
pytest
```

## Estrutura

```
Chaos/
├── src/chaos/        # Código-fonte do pacote
├── tests/            # Testes
├── pyproject.toml    # Configuração do projeto e dependências
└── README.md
```
