# Introduction

Welcome to my GitHub repository for the Triple DES (TDES) Encryption Class! This project is designed to implement the Triple Data Encryption Standard (3DES) algorithm, supporting Electronic Codebook (ECB), Cipher Block Chaining (CBC), and Output Feedback (OFB) modes.

## Project Overview

This repository contains a Python implementation of Triple DES encryption, a widely used symmetric-key cryptographic algorithm. The class is structured to allow easy encryption and decryption of messages using either an 8 byte or 24 byte key with multiple block cipher modes for enhanced security.

## Features

Modes of Operation: Supports **ECB**, **CBC**, and **OFB**

**Key Handling**: Accepts 8-byte (single key) or 24-byte (three keys) encryption keys

**Permutation Tables**: Implements initial/final permutations, S-Box substitution, and key scheduling

**Padding Support**: Ensures correct block alignment with PKCS#5 padding

**Unit Testing**: Includes comprehensive unit tests to validate functionality

## How It Works

**Key Expansion**: The provided key(s) are expanded into subkeys used in encryption rounds.

**Initial Permutation**: The plaintext undergoes an initial bitwise permutation.

**Feistel Rounds**: The data is processed through 16 different rounds of encryption with S-Box substitutions and XOR operations.

**Final Permutation**: The encrypted data undergoes a final permutation to produce the ciphertext.

**Decryption**: The same process is reversed using inverse subkeys.

Example Usage

``` python
# Initialize TDES with a 24-byte key in ECB mode
tdes = TDES(b'1234567890abcdef12345678', mode='ECB')

# Encrypt a message
plaintext = b'Hello, World!  '
ciphertext = tdes.encrypt(plaintext)
print(f'Ciphertext: {ciphertext.hex()}')

# Decrypt the message
decrypted = tdes.decrypt(ciphertext)
print(f'Decrypted: {decrypted.decode()}')
```
Author

This implementation was developed by Jordan Lassers as part of a cryptography course project.
