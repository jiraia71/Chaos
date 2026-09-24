"""Funções principais do projeto Chaos."""


def saudacao(nome: str = "mundo") -> str:
    """Retorna uma mensagem de saudação.

    Args:
        nome: nome a saudar. Padrão é "mundo".

    Returns:
        A mensagem de saudação.
    """
    return f"Olá, {nome}!"
