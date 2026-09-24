"""Testes para chaos.core."""

from chaos.core import saudacao


def test_saudacao_padrao():
    assert saudacao() == "Olá, mundo!"


def test_saudacao_com_nome():
    assert saudacao("Chaos") == "Olá, Chaos!"
