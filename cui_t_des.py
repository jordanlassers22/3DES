class TDES:
    """ Implements the Triple DES algorithm with a 192-bit key and three block
        modes: ECB, CBC, and OFB. """       
    
    # 32-bit to 48-bit
    _EXPAND = [31,  0,  1,  2,  3,  4,  3,  4,
                5,  6,  7,  8,  7,  8,  9, 10,
               11, 12, 11, 12, 13, 14, 15, 16,
               15, 16, 17, 18, 19, 20, 19, 20,
               21, 22, 23, 24, 23, 24, 25, 26,
               27, 28, 27, 28, 29, 30, 31,  0]
    
    # 32-bit permutation after S-BOX substitution
    _SBOX_PERM = [15,  6, 19, 20, 28, 11, 27, 16,
                   0, 14, 22, 25,  4, 17, 30,  9,
                   1,  7, 23, 13, 31, 26,  2,  8,
                  18, 12, 29,  5, 21, 10,  3, 24]
    
    # Initial permutation on incoming block
    _INIT_PERMUTATION = [57, 49, 41, 33, 25, 17,  9, 1,
                         59, 51, 43, 35, 27, 19, 11, 3,
                         61, 53, 45, 37, 29, 21, 13, 5,
                         63, 55, 47, 39, 31, 23, 15, 7,
                         56, 48, 40, 32, 24, 16,  8, 0,
                         58, 50, 42, 34, 26, 18, 10, 2,
                         60, 52, 44, 36, 28, 20, 12, 4,
                         62, 54, 46, 38, 30, 22, 14, 6]
    
    # Inverse of _INITIAL_PERMUTATION
    _FINAL_PERMUTATION = [39,  7, 47, 15, 55, 23, 63, 31,
                          38,  6, 46, 14, 54, 22, 62, 30,
                          37,  5, 45, 13, 53, 21, 61, 29,
                          36,  4, 44, 12, 52, 20, 60, 28,
                          35,  3, 43, 11, 51, 19, 59, 27,
                          34,  2, 42, 10, 50, 18, 58, 26,
                          33,  1, 41,  9, 49, 17, 57, 25,
                          32,  0, 40,  8, 48, 16, 56, 24]
    
    _S_BOXES = [
        [[14,  4, 13,  1,  2, 15, 11,  8,  3, 10,  6, 12,  5,  9,  0,  7],
         [ 0, 15,  7,  4, 14,  2, 13,  1, 10,  6, 12, 11,  9,  5,  3,  8],
         [ 4,  1, 14,  8, 13,  6,  2, 11, 15, 12,  9,  7,  3, 10,  5,  0],
         [15, 12,  8,  2,  4,  9,  1,  7,  5, 11,  3, 14, 10,  0,  6, 13],
        ],
        [[15,  1,  8, 14,  6, 11,  3,  4,  9,  7,  2, 13, 12,  0,  5, 10],
         [ 3, 13,  4,  7, 15,  2,  8, 14, 12,  0,  1, 10,  6,  9, 11,  5],
         [ 0, 14,  7, 11, 10,  4, 13,  1,  5,  8, 12,  6,  9,  3,  2, 15],
         [13,  8, 10,  1,  3, 15,  4,  2, 11,  6,  7, 12,  0,  5, 14,  9],
        ],
        [[10,  0,  9, 14,  6,  3, 15,  5,  1, 13, 12,  7, 11,  4,  2,  8],
         [13,  7,  0,  9,  3,  4,  6, 10,  2,  8,  5, 14, 12, 11, 15,  1],
         [13,  6,  4,  9,  8, 15,  3,  0, 11,  1,  2, 12,  5, 10, 14,  7],
         [ 1, 10, 13,  0,  6,  9,  8,  7,  4, 15, 14,  3, 11,  5,  2, 12],
        ],
        [[ 7, 13, 14,  3,  0,  6,  9, 10,  1,  2,  8,  5, 11, 12,  4, 15],
         [13,  8, 11,  5,  6, 15,  0,  3,  4,  7,  2, 12,  1, 10, 14,  9],
         [10,  6,  9,  0, 12, 11,  7, 13, 15,  1,  3, 14,  5,  2,  8,  4],
         [ 3, 15,  0,  6, 10,  1, 13,  8,  9,  4,  5, 11, 12,  7,  2, 14],
        ],
        [[ 2, 12,  4,  1,  7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9],
         [14, 11,  2, 12,  4,  7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6],
         [ 4,  2,  1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14],
         [11,  8, 12,  7,  1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3],
        ],
        [[12,  1, 10, 15,  9,  2,  6,  8,  0, 13,  3,  4, 14,  7,  5, 11],
         [10, 15,  4,  2,  7, 12,  9,  5,  6,  1, 13, 14,  0, 11,  3,  8],
         [ 9, 14, 15,  5,  2,  8, 12,  3,  7,  0,  4, 10,  1, 13, 11,  6],
         [ 4,  3,  2, 12,  9,  5, 15, 10, 11, 14,  1,  7,  6,  0,  8, 13],
        ],
        [[ 4, 11,  2, 14, 15,  0,  8, 13,  3, 12,  9,  7,  5, 10,  6,  1],
         [13,  0, 11,  7,  4,  9,  1, 10, 14,  3,  5, 12,  2, 15,  8,  6],
         [ 1,  4, 11, 13, 12,  3,  7, 14, 10, 15,  6,  8,  0,  5,  9,  2],
         [ 6, 11, 13,  8,  1,  4, 10,  7,  9,  5,  0, 15, 14,  2,  3, 12],
        ],
        [[13,  2,  8,  4,  6, 15, 11,  1, 10,  9,  3, 14,  5,  0, 12,  7],
         [ 1, 15, 13,  8, 10,  3,  7,  4, 12,  5,  6, 11,  0, 14,  9,  2],
         [ 7, 11,  4,  1,  9, 12, 14,  2,  0,  6, 10, 13, 15,  3,  5,  8],
         [ 2,  1, 14,  7,  4, 10,  8, 13, 15, 12,  9,  0,  3,  5,  6, 11],
        ]
    ]
    
    
    # 64-bit to 56-bit permutation on the key
    _KEY_PERMUTATION1 = [56, 48, 40, 32, 24, 16,  8,  0, 
                         57, 49, 41, 33, 25, 17,  9,  1,
                         58, 50, 42, 34, 26, 18, 10,  2, 
                         59, 51, 43, 35, 62, 54, 46, 38, 
                         30, 22, 14,  6, 61, 53, 45, 37,
                         29, 21, 13,  5, 60, 52, 44, 36,
                         28, 20, 12,  4, 27, 19, 11,  3]
    
    # 56-bit to 48-bit permutation on the key
    _KEY_PERMUTATION2 = [13, 16, 10, 23,  0,  4,  2, 27,
                         14,  5, 20,  9, 22, 18, 11,  3, 
                         25,  7, 15,  6, 26, 19, 12,  1,
                         40, 51, 30, 36, 46, 54, 29, 39, 
                         50, 44, 32, 47, 43, 48, 38, 55, 
                         33, 52, 45, 41, 49, 35, 28, 31]
    
    # Matrix that determines the shift for each round of keys
    _KEY_SHIFT = [ 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]
    
    def __init__(self, key, mode = "ECB", iv=b'\x00\x00\x00\x00\x00\x00\x00\x00'):
       """ Class Constructor.
   
           Parameters:
             key (bytes):  64-bit key used for DES encryption
   
           Returns:
             nothing
       """
       if mode.upper() not in ("ECB", "CBC", "OFB"):
           raise ValueError(f"{mode} is not a valid option. Valid options are: ECB, CBC, OFB")
       else:
           self._key = key
           self._mode = mode.upper()
           self.iv = iv
           self._iv = TDES._bytes_to_bit_array(self.iv)   
           self._split_into_subkeys()
    def encrypt(self, plaintext):
        """ Encrypts plaintext data with TDES (Data Encryption Standard).
    
            Parameters:
              data (bytes): input data to be encrypted
    
            Returns:
              An encrypted byte string of equal length to the original data
        """
        self.reset()
        
        if self._mode == "ECB" or self._mode == "CBC":
            plaintext = self._add_padding(plaintext)
        subkeys1 = self._generate_subkeys(self.subkey1)
        subkeys2 = list(reversed(self._generate_subkeys(self.subkey2)))
        subkeys3 = self._generate_subkeys(self.subkey3)
        plaintext = self._bytes_to_bit_array(plaintext)
        ciphertext = []
        
        if self.mode == "ECB":
            for plaintextBlock in self._nsplit(plaintext, 64):
                ct_block = self.encrypt_block(plaintextBlock, subkeys1)
                ct_block = self.encrypt_block(ct_block, subkeys2)
                ct_block = self.encrypt_block(ct_block, subkeys3)
                ciphertext += ct_block
                 
        elif self.mode == "CBC":
            for plaintextBlock in self._nsplit(plaintext,64):
                plaintextBlock = self._xor(plaintextBlock, self._iv)
                ciphertextBlock = self.encrypt_block(plaintextBlock, subkeys1)
                ciphertextBlock = self.encrypt_block(ciphertextBlock, subkeys2)
                ciphertextBlock = self.encrypt_block(ciphertextBlock, subkeys3)
                self._iv = ciphertextBlock
                ciphertext += ciphertextBlock
                
        elif self.mode == "OFB":
            for plaintextBlock in self._nsplit(plaintext,64):
                self._iv = self.encrypt_block(self._iv, subkeys1)
                self._iv = self.encrypt_block(self._iv, subkeys2)
                self._iv = self.encrypt_block(self._iv, subkeys3)
                result = TDES._xor(plaintextBlock, self._iv)
                ciphertext += result
        
        else:
            raise ValueError("Mode must be ECB, CBC, or OFB")
        ciphertext = TDES._bit_array_to_bytes(ciphertext)
        
        return ciphertext
    
    def decrypt(self, ciphertext):
        """ Decrypts given cyphertext with TDES (Data Encryption Standard).
    
            Parameters:
              ciphertext (bytes): input data to be decrypted
    
            Returns:
              A decrypted byte string
        """
        self.reset()
        # Generate the 16 subkeys
        subkeys = self._generate_subkeys(self._key) #remove this line once done
        if self.mode == "ECB" or self.mode == "CBC":
            subkeys1 = list(reversed(self._generate_subkeys(self.subkey1)))
            subkeys2 = self._generate_subkeys(self.subkey2)
            subkeys3 = list(reversed(self._generate_subkeys(self.subkey3)))
            subkeys.reverse()
        else:
            subkeys1 = self._generate_subkeys(self.subkey1)
            subkeys2 = list(reversed(self._generate_subkeys(self.subkey2)))
            subkeys3 = self._generate_subkeys(self.subkey3)
            
        ciphertext = self._bytes_to_bit_array(ciphertext)
        result = []
        
        if self.mode == "ECB":
            for ct_block in self._nsplit(ciphertext, 64):
                pt_block = self.encrypt_block(ct_block, subkeys3)
                pt_block = self.encrypt_block(pt_block, subkeys2)
                pt_block = self.encrypt_block(pt_block, subkeys1)
                result += pt_block
        elif self.mode == "CBC":
            for ciphertextBlock in self._nsplit(ciphertext, 64):
                decryptedBlock = self.encrypt_block(ciphertextBlock, subkeys3)
                decryptedBlock = self.encrypt_block(decryptedBlock, subkeys2)
                decryptedBlock = self.encrypt_block(decryptedBlock, subkeys1)
                plaintext = TDES._xor(self._iv, decryptedBlock)
                self._iv = ciphertextBlock
                result += plaintext
                
        elif self.mode == "OFB":
            for ciphertextBlock in self._nsplit(ciphertext, 64):
                self._iv = self.encrypt_block(self._iv, subkeys1)
                self._iv = self.encrypt_block(self._iv, subkeys2)
                self._iv = self.encrypt_block(self._iv, subkeys3)
                plaintextBlock = TDES._xor(ciphertextBlock, self._iv)
                result += plaintextBlock
        else:
            raise ValueError("Mode must be ECB, CBC, or OFB")
        # Convert the bit array back to bytes
        result = self._bit_array_to_bytes(result)
        
        # Remove padding
        if self.mode == "ECB" or self.mode == "CBC":
            result = self._rem_padding(result)
        
        return result
    
    def _split_into_subkeys(self):
        
        if len(self._key) == 8:
            self.subkey1 = self._key[0:8]
            self.subkey2 = self._key[0:8]
            self.subkey3 = self._key[0:8]
            
        elif len(self._key) == 24:
            self.subkey1 = self._key[0:8]
            self.subkey2 = self._key[8:16]
            self.subkey3 = self._key[16:24]
        else:
            raise ValueError("key must be 8 or 24 bytes")
        
            
    
    def reset(self):
        """
        Resets the IV to the default or user entered value when declaring object to allow for another encryption or decryption
        """
        self._iv = TDES._bytes_to_bit_array(self.iv)
        return 
    
    @staticmethod
    def _permute(block, table):
        """
        Permutes the input block based on the provided permutation table.
    
        Parameters:
        block (list): The input list of bits to be permuted.
        table (list): The permutation table, which specifies the new order of the bits.
    
        Returns:
        list: A new list with bits permuted according to the table.
        """
        return [block[i] for i in table]
    
    
    
    @staticmethod
    def func_f(r, subkeyi):
        """ Function that assists with a section of TDES.
    
            Parameters:
              right side of ciphertext, subkey
    
            Returns:
              permuted and substituted cyphertext
        """
        tmp = TDES._permute(r,TDES._EXPAND)
        tmp = TDES._xor(tmp, subkeyi)
        tmp = TDES._substitute(tmp)
        tmp = TDES._permute(tmp, TDES._SBOX_PERM)
        return tmp
       
    @staticmethod
    def encrypt_block(block, subkeys): #both lists of 0s and 1s
        """ Encrypts a 64 bit portion of cyphertext.
    
            Parameters:
              block: binary list
              subkeys: 2d binary list of generated subkeys
    
            Returns:
              encrypted block of text
        """
        block = TDES._permute(block, TDES._INIT_PERMUTATION)
        left = block[:32]
        right = block[32:]
        for i in range(16):
            tmp = TDES.func_f(right, subkeys[i])
            tmp = TDES._xor(tmp, left)
            left = right
            right = tmp
        block = TDES._permute(right + left, TDES._FINAL_PERMUTATION)
        return block
         
        
    @staticmethod
    def _substitute(bit_array):
        """
        The S-BOX function receives a 48-bit input that is stored as a list of ones and zeroes. This 48-bit input
    is divided into eight 6-bit chunks. Each chunk is divided into a row and column number and used to
    lookup/substitute a 4-bit value from an S-BOX table. The first 6-bit chunk looks up a value in the first
    S-BOX table and produces the first 4 bits of the output. The second 6-bit chunk looks up a value in the
    second S-BOX table and produces the next 4 bits of the output. This pattern continues for all of the
    remaining chunks. Notice that we are substituting eight 4-bit chunks for what was originally eight 6-bit
    chunks. This means that our input will be 48 bits in size but the output will only be 32 bits long
        """
        if len(bit_array) != 48:
            raise ValueError("The bit array must be 48 bits long.")
        output = []
        tableNum = 0
        for i in range(0, len(bit_array), 6):
            
            chunk = bit_array[i:i+6]
            #turn int list into strings
            rowNum = int(f"{chunk[0]}{chunk[5]}",2)
            colNum = int(f"{chunk[1]}{chunk[2]}{chunk[3]}{chunk[4]}",2)
            
            sboxVal = TDES._S_BOXES[tableNum][rowNum][colNum] #look up value in sbox tables
            binString = bin(sboxVal)[2:].zfill(4) #keep leading 0s, remove prefix
            binList = [int(bit) for bit in binString]
            output += binList
            tableNum += 1
        return output
    
    @staticmethod
    def _generate_subkeys(key:bytes):
        """
        Generates 16, 48bit subkeys for TDES implementation
        Takes in a 64 bit key
        """
        subkeys = []
        keybits = TDES._bytes_to_bit_array(key)
        k0 = TDES._permute(keybits, TDES._KEY_PERMUTATION1) #Permutes from 64 bits to 58 bits
        right = k0[28:]
        left = k0[:28]
        for i in range(16): #creates 16 subkeys bit shifting by 1 to the left for each one
            left = TDES._lshift(left, TDES._KEY_SHIFT[i])
            right = TDES._lshift(right, TDES._KEY_SHIFT[i])
            ki = TDES._permute(left+right, TDES._KEY_PERMUTATION2)
            subkeys.append(ki)
        return subkeys
    
    @staticmethod
    def _add_padding(message):
        """
        Adds padding to ensure message is a multiple of 8
        """
        padding_length = 8-len(message)%8
        padding = chr(padding_length) * padding_length
        message += padding.encode("utf-8")
        return message
    
    @staticmethod
    def _rem_padding (message):
        """Removes padding from a message padded using the self._add_padding function."""
        padding_length = message[-1]  # Get the last byte which represents the padding length
        return message[:-padding_length]  # Slice off the padding
    
    @staticmethod
    def _bytes_to_bit_array(byte_string):
        """Convert a bytes object into a list of bits"""
        result = []
        for byte in byte_string:
            bits = [(byte >> i) & 1 for i in range(7, -1, -1)]
            result.extend(bits)
        return result  
    
    @staticmethod
    def _bit_array_to_bytes(bit_array):
        """Convert a list of bits into a bytes object"""
        if len(bit_array) % 8 != 0:
            raise ValueError("Bit array length must be a multiple of 8.")
        
        byte_array = bytearray()
        for i in range(0, len(bit_array), 8):
            byte = sum([bit << (7 - j) for j, bit in enumerate(bit_array[i:i+8])])
            byte_array.append(byte)
        
        return bytes(byte_array)
    
    @staticmethod
    def _nsplit(data, split_size=64):
        """Splits data into blocks of `split_size` bytes. 
           If the data length is not divisible by `split_size`, the last block will be smaller. """
        result = []
        for i in range(0, len(data), split_size):
            result.append(data[i:i + split_size])
        return result
    
    @staticmethod
    def _hex_print(block):
        s = [str(i) for i in block]
        b = int("".join(s),2)
        print(hex(b)[2:].zfill(16).upper())
        return
    
    @staticmethod
    def _lshift(sequence: list, n: int):
        """
        Left shifts the elements in the sequence by n positions with wrapping.
    
        Parameters:
        sequence (list): The input list of bits or elements to be shifted.
        n (int): The number of positions to shift.
    
        Returns:
        list: A new list where the elements are shifted to the left by n positions.
        """
        length = len(sequence)
        n = n % length  # Ensure n is within the length of the sequence
        return sequence[n:] + sequence[:n]
    
    @staticmethod
    def _xor(x, y):
        """
        XORs two lists element by element.
    
        Parameters:
        x (list): The first input list of bits.
        y (list): The second input list of bits.
    
        Returns:
        list: A new list where each element is the result of XORing the corresponding elements of x and y.
        
        """
        return [a ^ b for a, b in zip(x, y)]
    
    # Getter for key
    @property
    def key(self):
        return self._key
    
    # Setter for key
    @key.setter
    def key(self, new_key):
        if len(new_key) == 8 or len(new_key) == 24:  # Ensure the key is 8 bytes or 24 bytes long (64 bits)
            self._key = new_key
            self._split_into_subkeys()
        else:
            raise ValueError("Key must be 8 bytes or 24 bytes long.")
        
    @property
    def mode(self):
        return self._mode
    
    # Setter for mode
    @mode.setter
    def mode(self, new_mode):
        if new_mode.upper() not in ("ECB", "CBC", "OFB"):
            raise ValueError("Mode must be ECB, CBC, or OFB")
        self._mode = new_mode.upper()
    
def run_unit_tests():
    """ Runs unit tests for each function in this module. Prints 'ALL UNIT
        TESTS PASSED' if all of the unit tests were successful. Raises an
        AssertionError if any single unit test fails. """
    try:
        # self._add_padding
        test1result = TDES._add_padding(b'CSC428')
        assert test1result == b'CSC428\x02\x02', "Unit test #1 failed: self.self._add_padding(b'CSC428')"
        
        test2result = TDES._add_padding(b'TALLMAN')
        assert test2result == b'TALLMAN\x01', "Unit test #2 failed: self.self._add_padding(b'TALLMAN')"
        
        test3result = TDES._add_padding(b'JTALLMAN')
        assert test3result == b'JTALLMAN\x08\x08\x08\x08\x08\x08\x08\x08', "Unit test #3 failed: self.self._add_padding(b'JTALLMAN')"

        # _rem_padding
        test1result = TDES._rem_padding(b'CSC428\x02\x02')
        assert test1result == b'CSC428', "Unit test #1 failed: _rem_padding(b'CSC428')"
       
        test2result = TDES._rem_padding(b'TALLMAN\x01')
        assert test2result == b'TALLMAN', "Unit test #2 failed: _rem_padding(b'TALLMAN')"
        
        test3result = TDES._rem_padding(b'JTALLMAN\x08\x08\x08\x08\x08\x08\x08\x08')
        assert test3result == b'JTALLMAN', "Unit test #3 failed: _rem_padding(b'JTALLMAN')"
        
        # _bytes_to_bit_array
        test1result = TDES._bytes_to_bit_array(b'\x00')
        assert test1result == [0,0,0,0,0,0,0,0], "Unit test #1 failed: _bytes_to_bit_array(b'\x00')"
        
        test2result = TDES._bytes_to_bit_array(b'\xA5')
        assert test2result == [1,0,1,0,0,1,0,1], "Unit test #2 failed: _bytes_to_bit_array(b'\xA5')"
       
        test3result = TDES._bytes_to_bit_array(b'\xFF')
        assert test3result == [1,1,1,1,1,1,1,1], "Unit test #3 failed: _bytes_to_bit_array(b'\xFF')"
        
        
        # _bit_array_to_bytes
        test1result = TDES._bit_array_to_bytes([0,0,0,0,0,0,0,0])
        assert test1result == b'\x00', "Unit test #1 failed: _bit_array_to_bytes([0,0,0,0,0,0,0,0])"
        
        test2result = TDES._bit_array_to_bytes([1,0,1,0,0,1,0,1])
        assert test2result == b'\xA5', "Unit test #2 failed: _bit_array_to_bytes([1,0,1,0,0,1,0,1])"
    
        test3result =TDES. _bit_array_to_bytes([1,1,1,1,1,1,1,1])
        assert test3result == b'\xFF', "Unit test #3 failed: _bit_array_to_bytes([1,1,1,1,1,1,1,1])"
        
        #_nsplit
        test1result = TDES._nsplit(b'1111222233334444', 4)
        assert test1result == [b'1111', b'2222', b'3333', b'4444'], "Unit test #1 failed: _nsplit(b'1111222233334444', 4)"
        
        test2result = TDES._nsplit(b'ABCDEFGHIJKLMN', 3)
        assert test2result == [b'ABC', b'DEF', b'GHI', b'JKL',b'MN'], "Unit test #2 failed: _nsplit(b'ABCDEFGHIJKLMN', 3)"
        
        test3result = TDES._nsplit(b'THE CODE BOOK BY SINGH', 5)
        assert test3result == [b'THE C', b'ODE B', b'OOK B', b'Y SIN',b'GH'], "Unit test #3 failed: _nsplit(b'THE CODE BOOK BY SINGH', 5)"
        
        
        
        #Permutation Function Testing
        initPermTestBlock = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
 '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 
 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
 '!', '$']
        initPermTestResult = TDES._permute(initPermTestBlock,TDES._INIT_PERMUTATION )
        assert initPermTestResult == ['v', 'n', 'f', '8', 'Z', 'R', 'J', 'B', 'x', 'p', 'h', '0', '2', 'T', 'L', 'D', 'z', 'r', 'j', 'b', '4', 'V', 'N', 'F', '$', 't', 'l', 'd', '6', 'X', 'P', 'H', 'u', 'm', 'e', '7', 'Y', 'Q', 'I', 'A', 'w', 'o', 'g', '9', '1', 'S', 'K', 'C', 'y', 'q', 'i', 'a', '3', 'U', 'M', 'E', '!', 's', 'k', 'c', '5', 'W', 'O', 'G'], "Init Permute Test failed"
        
        
        finiPermTestBlock = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
 '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 
 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 
 '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 
 '!', '?', '*', ':', ')']
        finalPermTestResult = TDES._permute(finiPermTestBlock,TDES._FINAL_PERMUTATION )
        assert finalPermTestResult == ['d', 'H', 'l', 'P', '7', 'X', ')', '5', 'c', 'G', 'k', 'O', '6', 'W', ':', '4', 'b', 'F', 'j', 'N', '5', 'V', '*', '3', 'a', 'E', 'i', 'M', '4', 'U', '?', '2', '9', 'D', 'h', 'L', '3', 'T', '!', '1', '8', 'C', 'g', 'K', '2', 'S', '0', '0', '7', 'B', 'f', 'J', '1', 'R', '9', 'Z', '6', 'A', 'e', 'I', '0', 'Q', '8', 'Y'], "Final Permute Test failed"
        
        
        expandTestBlock = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
 '0', '1', '2', '3', '4', '!']
        expandTestResult = TDES._permute(expandTestBlock, TDES._EXPAND )
        assert expandTestResult == ['!', 'A', 'B', 'C', 'D', 'E', 'D', 'E', 'F', 'G', 'H', 'I', 'H', 'I', 'J', 'K', 'L', 'M', 'L', 'M', 'N', 'O', 'P', 'Q', 'P', 'Q', 'R', 'S', 'T', 'U', 'T', 'U', 'V', 'W', 'X', 'Y', 'X', 'Y', 'Z', '0', '1', '2', '1', '2', '3', '4', '!', 'A'], "Expand Permute Test failed"
        
        
        sBoxPermTestBlock = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
 '!', '1', '2', '3', '4', '5']
        sBoxPermTestResult = TDES._permute(sBoxPermTestBlock, TDES._SBOX_PERM )
        assert sBoxPermTestResult == ['P', 'G', 'T', 'U', '2', 'L', '1', 'Q', 'A', 'O', 'W', 'Z', 'E', 'R', '4', 'J', 'B', 'H', 'X', 'N', '5', '!', 'C', 'I', 'S', 'M', '3', 'F', 'V', 'K', 'D', 'Y'], "Sbox Permute Test failed"
        
        #L shift test
        lShiftTestResult1 = TDES._lshift([1, 0, 1, 1, 0, 0, 1], 3)
        assert lShiftTestResult1 == [1,0,0,1,1,0,1], "L shift Test 1 failed: _lshift([1, 0, 1, 1, 0, 0, 1], 3)"
        
        lShiftTestResult2 = TDES._lshift([1, 1, 1, 1, 0, 0], 2)
        assert lShiftTestResult2 == [1, 1, 0, 0, 1, 1 ], "L shift Test 2 failed: _lshift([1, 1, 1, 1, 0, 0], 2)"
        
        # XOR Test
        x = [1, 0, 1, 1, 0, 0, 1]
        y = [0, 1, 1, 0, 0, 1, 1]
        xORTestResult = TDES._xor(x, y)
        assert xORTestResult == [1, 1, 0, 1, 0, 1, 0], "XOR Test failed:  xORTestResult = _xor(x, y)"
        
        
        #Substitution test
        subTestBitArray = [1,0,1,1,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,0,0,1,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,1,1,0,0,1,0,0]
        subTestCorrectAnswer = [0,0,1,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,0,1,0,0]
        substituteTestResult = TDES._substitute(subTestBitArray)
        assert substituteTestResult == subTestCorrectAnswer, "Substitute test failed: _substitute(subTestBitArray)"
        
        #SubKey test
        subkey_input = b"\xEF\x00\xEF\x00\xFF\x80\xFF\x80"
        subkey_answer = [ [0,1,1,0,1,1,1,1,1,0,1,0,1,1,0,0,0,0,0,1,1,0,1,1,
                   1,0,1,1,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,0],
                  [1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,1,1,1,1,0,1,1,0,1,
                   0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1],
                  [1,0,0,1,0,0,0,1,0,1,0,1,0,0,1,1,1,1,1,0,1,1,0,1,
                   0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,1,1,0,1,1,0,1,0,1],
                  [1,0,0,1,0,0,0,1,0,1,0,1,1,0,1,1,1,1,1,0,0,1,0,1,
                   0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,1,0,1,1,0,1,0,1],
                  [1,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,1,1,1,0,0,1,0,1,
                   0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,1,0,0,1,1,1,0,1],
                  [1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,0,0,1,0,1,
                   0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0,1],
                  [1,1,0,1,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,0,0,1,0,1,
                   0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,1,1,0,1,0,1,1,0,1],
                  [1,1,0,1,0,0,0,1,1,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,
                   0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,1,0,1,1,0,1],
                  [1,1,1,0,1,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,1,0,
                   0,0,1,1,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,1,0],
                  [1,1,1,0,1,1,1,0,1,0,1,0,1,1,1,0,0,0,0,1,1,0,1,0,
                   1,0,1,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,1,0],
                  [0,1,1,0,1,1,1,0,1,0,1,1,1,1,1,0,0,0,0,1,1,0,1,0,
                   1,0,1,0,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,0,0,0,1,0],
                  [0,1,1,0,1,1,1,0,1,0,1,1,1,1,0,0,0,1,0,1,1,0,1,0,
                   1,0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0],
                  [0,1,1,0,1,1,1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,
                   1,0,0,1,1,1,0,0,1,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0],
                  [0,1,1,0,1,1,1,0,1,1,1,0,1,1,0,1,0,0,0,1,1,0,1,0,
                   1,0,0,1,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0],
                  [0,1,1,0,1,1,1,0,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,1,
                   1,0,1,1,1,0,0,0,1,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0],
                  [1,0,0,1,1,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,
                   0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,1,1,1,0,1] ]

        generateSubkeysResult = TDES._generate_subkeys(subkey_input)
        assert generateSubkeysResult == subkey_answer, "Generate subkey test failed: _generate_subkeys(subkey_input)"
        
        #Func F test
        rightTest = [1,0,1,1,1,1,1,0,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,1,0,1,0,1,1]
        subkeyTest = [1,0,1,0,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,1,0,1,1,0,0,0,0]
        funcFExpectedResult = [1,1,1,1,0,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,0]
        funcFResults = TDES.func_f(rightTest, subkeyTest)
        assert funcFResults == funcFExpectedResult, "func_f test failed: func_f(rightTest, subkeyTest)"
        
        # print("Hex Print output:")
        # _hex_print([1,1,1,1,0,1,0,1,0,0,0,0,1,0,1,0
        #            ,1,0,0,1,0,1,1,0,1,1,0,1,1,0,1,1])
        
        #ECB Encrypt Test 1
        testTDES = TDES(b"\xaa\xbb\xcc\xdd\xee\xff\x00\x11", "ECB")
        assert testTDES.encrypt(b"I can't believe this worked!").hex() == "359513b15ee9e9d37cb4bbfad05bc8645588a626d18f83e07c04cb1d798d2e2b", "ECB encryption test 1 failed: encrypt(plaintext, key)"
        
        #ECB Encrypt Test 2
        testTDES = TDES(b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9', "ECB") 
        assert testTDES.encrypt(b"I can't believe this worked!").hex() == "df27aacac35420845731575f2e826b6abbdf3787b472defe963c564b8c5f96a1", "ECB encryption test 2 failed: encrypt(plaintext, key)"
        
        #ECB Decrypt Test 1
        #testTDES = TDES(b"\xaa\xbb\xcc\xdd\xee\xff\x00\x11", "ECB")
        #assert testTDES.decrypt(b'5\x95\x13\xb1^\xe9\xe9\xd3|\xb4\xbb\xfa\xd0[\xc8dU\x88\xa6&\xd1\x8f\x83\xe0|\x04\xcb\x1dy\x8d.+') == b"I can't believe this worked!", "ECB decryption test 1 failed: decrypt(ciphertext, key)"
        
        #ECB Decrypt Test 2
        testTDES = TDES(b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9', "ECB") 
        assert testTDES.decrypt(b"\xdf'\xaa\xca\xc3T \x84W1W_.\x82kj\xbb\xdf7\x87\xb4r\xde\xfe\x96<VK\x8c_\x96\xa1") == b"I can't believe this worked!", "ECB decryption test 2 failed: decrypt(ciphertext, key)"
        
        #CBC Encrypt Test 1
        testTDES = TDES(b"\xaa\xbb\xcc\xdd\xee\xff\x00\x11", "CBC", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.encrypt(b"I can't believe this worked!").hex() == "359513b15ee9e9d3d778feaaa4dfdf7c637d0f500e9494cabbfd5d2f10f2a4f1", "CBC encryption test 1 failed: encrypt(plaintext, key)"
        
        #CBC Encrypt Test 2
        testTDES = TDES(b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9', "CBC", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.encrypt(b"I can't believe this worked!").hex() == "df27aacac3542084cad3cfbe37b3b7815faf76975cf47d60f0c6e1f696907656", "CBC encryption test 2 failed: encrypt(plaintext, key)"
        
        #CBC Decrypt Test 1
        testTDES = TDES(b"\xaa\xbb\xcc\xdd\xee\xff\x00\x11", "CBC", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.decrypt(b'5\x95\x13\xb1^\xe9\xe9\xd3\xd7x\xfe\xaa\xa4\xdf\xdf|c}\x0fP\x0e\x94\x94\xca\xbb\xfd]/\x10\xf2\xa4\xf1') == b"I can't believe this worked!", "CBC decryption test 1 failed: decrypt(ciphertext, key)"

        #CBC Decrypt Test 2
        testTDES = TDES(b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9', "CBC", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.decrypt(b"\xdf'\xaa\xca\xc3T \x84\xca\xd3\xcf\xbe7\xb3\xb7\x81_\xafv\x97\\\xf4}`\xf0\xc6\xe1\xf6\x96\x90vV") == b"I can't believe this worked!", "CBC decryption test 2 failed: decrypt(ciphertext, key)"
        
        #OFB Encrypt Test 1
        testTDES = TDES(b"\xaa\xbb\xcc\xdd\xee\xff\x00\x11", "OFB", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.encrypt(b"I can't believe this worked!").hex() == "535d3ebba8560d00eac831f491af59e1521cf1387c946b88e852bf4b", "OFB encryption test 1 failed: encrypt(plaintext, key)"
        
        #OFB Encrypt Test 2
        testTDES = TDES(b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9', "OFB", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.encrypt(b"I can't believe this worked!").hex() == "2d1c633e3744231d28fa2de290453ff71556213a3251a5bdcd514abb", "OFB encryption test 2 failed: encrypt(plaintext, key)"
        
        #OFB Decrypt Test 1
        testTDES = TDES(b"\xaa\xbb\xcc\xdd\xee\xff\x00\x11", "OFB", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.decrypt(b'S]>\xbb\xa8V\r\x00\xea\xc81\xf4\x91\xafY\xe1R\x1c\xf18|\x94k\x88\xe8R\xbfK') == b"I can't believe this worked!", "OFB decryption test 1 failed: decrypt(plaintext, key)"

        #OFB Decrypt Test 2
        testTDES = TDES(b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9', "OFB", b'\x00\x00\x00\x00\x00\x00\x00\x00')
        assert testTDES.decrypt(b'-\x1cc>7D#\x1d(\xfa-\xe2\x90E?\xf7\x15V!:2Q\xa5\xbd\xcdQJ\xbb') == b"I can't believe this worked!", "OFB decryption test 2 failed: decrypt(plaintext, key)"

        print("ALL UNIT TESTS PASSED")
    except AssertionError as e:
        print(str(e))
    
if __name__ == "__main__":
    
    run_unit_tests()
    
    # ECB - Electronic Code Book
    secret_key1 = b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9'
    ciphertext1 = b'.<\xef\xec\x03\x9d\xc6>\x03U\xde\xdc\xfe,f\xbe\xef\xe3\x15\x13\x8e+m\xb1D,\x10\x89\xf5g\xcaCh\x88.\xd3\xe3\xcf\xfb\xd1\x87\xc22U\xe4\x07\x02\x17\xe6)!\x06\x8c\xeeu\xc1\xc7\xa5!Yb\xe3!\xe7\x02\xeb\xd2h\x97/\x8aUf7\x12i\x1ez\x07\x82\xf2M\xcf\xf0&\xc9O\xf49:0\xfdv\x0cT.D\xddc\xfe\xd8t\x8b\xf4\xa1oq0G}\r\xb4b\x1c\x9a\xad\x98\xfe\xad\x8a\x1f)\x88g\x9d\xb4\xed\x01H\x05\x9c^\xd5\x84G\xb6\xa6N\xbeo\xdd\xa9:\xf3\x9e\x16\x9b\xd5S\xdc6\xed"\x08\x8eK3\x82\xe2\x04\xd8*\x96X$\xaf\xdc\xa7>\x9f\'\xef\x88\xba\xca\x9b\x9cMn+\x91~W\xde\xcfx\x811^<\x8aS,\x01f\t\x87j4\xb1?\xef\x99\xa0\x98\x0cS\xe5\xb7\xcb\x1e\xb3\x0b\x10\x86L\xa2\x02\x83\x830\xb1\x17\xd6W\xbeB\x0f\x8e9 "Z\xad\xa3\xaa&\x98^\x1d\x01mn\xa1<\xb0\xab\xf8[\x06\xd6\x17\x94\xd7\xbdk\x98\xb8\xfc\xa56\t\x10\xb0\xba\x9f3\xe5\xe9\x8a\x98,-\x08`\xa8Q}\xcd\x99\x99\xea\xaa\xb8\xfa\x82\xda\xca\x9c_\x9ag\xee\xdf\x03\xcc\x19d$S\xc2\xbcsj\x99h\xcd2\xb4JpA\xa6\x1b\xce\xa1\xd2\xe6\xc1\x82\xd9ux{\xa5\xa8\xae\xdb\xc7\xf2\xa1\x027\xbbg\xf7\xfc\xd7\x08\x9c+Ks\xf9\x16\xf9/\xd2\x94"\x99\xee\xd52<\x9a\xec\x7f%\xef\xce\xbbu\x05\x88\x9f\x087\xe8\xd3(~6\na-\xa0\xa4{\x83\xaeza^3\xd8m`\x9bl\xc3\xc1\xbdZ\x9e\xf5\xab0\xef*\xdd:\xb2|u\n\xf1K\xac\xdc>c\x1c\x9e:\rfnP-\x08\xdd\x96\x9e\x7f^ \x90\xb7\x16+\x96\x1f\xca\xc4\x98\xa8h\xf0r\n\x0brn\xf5Qa\x06\xb23\xa0i]\xd8\xed\xbd\xe0\xda\xb3Zx\x97fh\x83\xd1\'d\xe4\x99\xa6\x94\xb6L*\x0c;\x8b\xadj]I\x81\x1e1\x99\x15\xff\xd0\xc4\xa9\xac\x08M\xe1)\x06\xebV\x89q\x84\xa5\x88F\x12\x95\xa0e7\xb3\xbc|L\xf7\x01\xf3\x9f\xa02z\xcf\xca\xa9\xbd\xccQ\xb6\xba\x95\xb3\xfb\x8dO!W'
    tDesObj = TDES(secret_key1)
    print(tDesObj.decrypt(ciphertext1).decode("utf-8"))
    print("\n")
    
    # CBC - Cipher Block Chaining
    secret_key2 = b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9'
    initvector2 = b'\x88\x99\xaa\xbb\xcc\xdd\xee\xff'
    ciphertext2 = b'\x9e\x8f\xe5G~\xd3\xeb\x9c\xd5\xe8\xb8\x1a\x17\xe5\xa1H\xebVX\xbe\xcf}W\xc8\xb5\xa6\xfb\xbd\x1fG\x8b\x13.\x1c*\x1d\t\xab\x1c\x8d\xf3#\xf8\x86\xf4\xbad\x18\x9fR\xad\xc7\xdd\xe5\xde\xf7\xd9\xca_\x1eY\xe4\x06-s\x86\x01\x0e\xdad\x15Z\x08\x92=6\x9f\xbc\x85\x8b\xfe\xd1\x9a/W 0\x16\xa1\x07\x04\x9dY\x8c\x85Z\xeb\x11r\xa3\xff\xb1\xbbd\xd9\xb8\xa7\x1d\x88\xdc;\xfb\x8d\xde\x04\x17\xb8\xde\xd5\xce*a\x93bQ`~\x17R\xec\xb928\x9f\x9d\xe4\xaa\x08.\xe9\xfdV\xf0\xa0\xa7`\xa1\xff^\x7f\xed,/El\xab%\x85\x91\xc3\xe2s\\\x96H\xc3\xceh\x8c\x91\x87\n\xa2u\x1b.@B!\xfbg71~)P\xa5\x9f\xd3R\xdb\xa1\xb33\xe7\xbc\xbc\xf2\x9d\xbf\x08\x99\xa3\xaf\xba\x16\x0eEh\xeb(/\x8cF`K\xe7\xaav\xe1\x9d\xc1\x98\xff\xf5r)e\x06\xa1\xa1\xef\x83W\xf9\xf3\t\xc4\x14p_:\xcc]\xca,\x11y\xdbf\xae\x9f\x7f\xb1Z\x95\\\xe3\xd4\x83A\xce\xc9\xe5\x97\x82pHv\x8d\x81\xcc\x9el\xa9}\x85A\xf0\xaf\x01B2\xd7\xc9\xe7Eff\xb5GY\x85\xd0\x01\xd2\xc9\x88N\x93\x86\x9a\xab\\T\xc2\xb9\xf6\\\\H\xed\xb0\xb5\xf8\x1azW\xc2K\x12\xfe\xc1\xd4\x86\x8c\xc9\xa6IU\x0c\xc4\xfa:\xd9c\x0c\xf4\x06\xf2Kg\x82D\x15\xc1_bG\x97\xbe\x9ed\xe8\x88\x87IS\xc1\xa6\x12\xb5\xa5\xe1.\xc5T\x191\n\xca\xa5\xfb\xe0\x92L`MG\xbd\xfa9\r\x7f\xb7U$^<\x1fx\xd3\xf6\xaa\xfaK\xbe\x81\xa6\x0e\xe5\x15#L\xfa\x1c\xc0\xa3\x08r\xb0~\xf6Jh\xf5\xf8\xfdx 3\xd2\xed\x87\xc2\x1c\x01\xea\xe9v]Vcl`\x8a\xb5\xcf\x08\xd0z\x00\xa3/\xc6\xa0\xaei*\x8b\xae\xd3\xb7\x87\x98B\xf2\x15\x08/N$\x80\x80\x0f6\x1c\r\x05\x12#:E\xd8\xf0vG\rv\xf8X\'\xc3\xa3\xc4\xda\xab,\xac\xecT*\x83\xc9\xd2tz(t\xe0\xa3d\xeb*\xc0D\xea\xbe\xbc\x18\xa4\xc7\xa0%a\x01G\x13\xa4\xb3\x12\xff!\xc2#\xee\xb7\xd4M\xfa\x9fTG2v\x14T\xd7\r\x1fJ\x90&\x80\x0fblPtG3`\xfe|~Q\x07h\r[$\xcdT:/\x94\xca\xce\x84\xd5\xa5@\xa0\xa6\x83\xba\x1b\xdf\x13~:hr\xa5\xad\x8b\x14\xfe\x1ae\xf3\x91\xe3\xd1\xe3\x99\xc0\xe9\xb7\xab\xac\x7f3i\x18\x8a\x1a"\xee\xc0Y<G\xea\xbf6\xba}2\x9a\xca\xf2OTM\xdeUQ\x13\xf3+\xdc\xcdj8R\xbb\x8a\xad\x1a\x16F\xf4`!\xc2\xeb\x12o\xf2\xb5\xefU\xc9\x01\x85l\xfe>\x15uF\x84\x1f\r\xfc\x8a{u\tLZ\xcf\xbd\x17\xb3X\xca\xc0\xa2\xaf\xb1\xb0\xb0>\x04\xdaK\x14\xb3\x1f\xd1\xabU\xb5z\xe4Z=\xa8\xf4\xc7\x9f[\xaf\xa2w\x9e#}\x9cO\xfa\xa8\xdb%\x8am\xfb1\xc2x\xa2\x9f1G\x02b\xa4c\xc7r\xa7yr\xe1-\xbc\xdfq\xcc=\x028\xc4\r\xf1x\x8bzc\xff`h\xdf(\x8eRS\xde\x8e\xa5\xc2<n\xcf\xfe\xe0\x8e\xb1\xd8\xc3N\xe9]\x1a\xe6\xaf\xf3DB\xc2\xa4:\xad>\xbd\xc2\xb1\x89\xe5\x98\xdb]\x14\x91\xf9\xdf\xcf\x11\x1f\xc5p\rr\xe0\xf8\x19S\x9f\xee\x08\r\x02\xafF/\xc7\x9a5UPb\x1a\x9a\x13Yv\xe3\xf9\xb3\x1d\xb5\x93\xe5\\HH\xc7\xe2\xd6}\x87\xc57g\x17\xa7\xbe\x96\x80"9u\x9a\x9e\xdc\xd7\xb8\xff\x96\xab\xa9\xd0\xee\x11CzIc\x16\xd9\x98X\x14\xe7\xcb\x89k\r4\x05@\x87\xc9\xf78\xb8\x9c\xba\xad+v\xee.\xacON\x94\xc8\xf3\xb7)\xfe\xcb\x97p\r\x94\x95\xd7\n\xdf\x144\x1c>\x0f\xfd\x08]\xa2\xb4\x05\n\xb5\xea\x0b\xeb\x11\x9e\x04A\x1c;\xb6\x00\xc9\xf3\x8d9K\x12\x7fE\x98\xeb\x10MM \x95)5\x9b\x01r\x9a\x96\x17m\xb4\xf7\xadC\xa5f\xf8y\n\x0c.\xd3\xf2\x06\x14Q_\xdf\x11\x80y\x0c\n\xccuu\x87\xeb~\x9f\xc4\xde\x83\x96a\xdc\xbc\xf6N\r?,\xf7+k\xde\xc7\x04\x1fJ\xc9&\x91z\xc8\xd5M{\xa6\xad\xc9\x9d\xcfF\x89^\xe3\xb8S\xabQ\xd1\xf4\x17\xc6\xd3\x99\xf4\xe8\xab1s\xc2A\x96\xdd\xd8i\x9f\xf6>Hu\x8b\xf7 "\x8c\xdf\xc7$U;\t2\x1ce\x1c\xf9l\x1d+\xe5\x05\xd3d\x98?\xc4\x14=o\xbb \xd9\x01\xb3\xbb\xa2\xc7o\xc1<f|\x8e\xb1^\x95\xf6\x98\x03dy\x1da\x94uH\x94\xeb\x039d\xd9\x107\x13\x16\x80\x7f\xd3X\xd4\x95\x070ye\xc75\x8e*de\x06\x05U+\x1d\xdf\xf3\xf1\xb5\x06K\xf4\x7f"{ \xdd5|\xd2G*\xd2\x06\xc8\x10\xac\xed\xe2\x1f=\x06'
    tDesObj = TDES(secret_key2, "CBC", initvector2)
    print(tDesObj.decrypt(ciphertext2).decode("utf-8"))
    print("\n")
    
    # OFB - Output Feedback
    secret_key3 = b'\xde\xad\xbe\xef\x7a\xc0\xba\xbe\xca\xfe\xf0\x0d\xca\x7d\x00\xd1\x23\x45\x67\x89\x0a\xbc\xde\xf9'
    initvector3 = b'\x88\x99\xaa\xbb\xcc\xdd\xee\xff'
    ciphertext3 = b'\xaa\xa6\\6\x9f\xbb)\x1b\xb2\x14\xfb\xf2\x88\xebQ?\xa3\xaf\x8b|l\xdaG>\xd1\x84\x83l(\x18\x9a\x1cG\x11t\x9b\xf1\x9b<\r)SD \xd0\x01\x05\xc5\xdb\xed\x8e\x1a)\xbd\x8c\xfcG\xbfB\xf2\xc8%R\x04\x8a\x92{\xc1VPV\xea\x98\x1e\x0e\xfaA\xcb\xb6!\x96]wA\x98\xbe*\x86\xb3\x1a\xdc$\x85x\xe1HM\xf5Z\xcd\x10\xfc\xa6d\x17\x0f#\xb7\xdcu\xbe\x04B46\x1a69c\x10&\x0f\xb0,\x85\xb8\x04\xa7\x9f\x95\xe3\x06\xf3\x17L\xd7\xff\x8c\xee\xf42\xb7\xa9\xfcJuB7\x83\x83\x8b\x89Eq8\x80\xb4\x7fQ\xa8\xd089O\xf6\xc81j\xc3u\xb3[\xcb\xe5,P\xebqT\xe5\xf9\xc1\x98\xd3\x1b\x9c\xb6\'v&\x14\x86-\x1c^\xad\x7f\xd2\xcb\xb9t(gh\xb9\xc5\x9c$\xdc@\x86\xa1\xc2\x98v\x9a\r\x05\xef\x99\x89\xe3\xcb\xe6\x1a\x89t\x04;\xb5@@E\xc2\x15\x9a\x9f6\xaf8\x9c\x04E\x08Qf\xb7\xbeD\xdc\xf8\xbag*?\xd5Zl\xa4\xabW\xcc[S\xe9"\x03c\xc7\xb0\r\xe0\x17+(.\x85?\x06\xbak\xdc\x8b\xba6A-\xf7c&\xc2\xbets\xb0\x0e\xe7h\xe8\xd6\r\x1a\x8d\xd7\xad8\x8b\xf4h\xdc\x15\xe3\x04\xf6\xbc\xe2\xae\x04\n\t\r\x83(k4\x0e\x15\x90\xeb\xc9z\xfc>\x14H\xf0S5P\xd3\x86\xacpn\x1f\xcc\xeej\x1b_\xeb\x1f\xcc\xeb\x01\xad\x16G\xab2w\xc0\xa7\xcb\xd6\xbd\xad9\xd5\x1f\x1cu\x9c6l\x9aQ\x8f\x8cu#\xb829s_\x05h\x90\x847J\xa4\x87`\xc9\xb0\x8d\xf8HGRw\xe6e\xbe\xb5\xfc\xa2\xfc\x1c\xc6"\xb7\x11%\x8c#\xe3\xeat\xf9\xf0\x1f\x08>o\x11\xe4X)\x8e\x08\x1e:\xc9qP\xc0p\xe5t\x98\xa7l\x91\x0c\x94\xe6c\x98E\xac\xfb\x13p\xd5\xc3/et- \xa6\xa2\x9d\x04\x18\xe0\x0b\xea\xe1\x96\xec\xb6^\x8d\x99\x80\xc2\x93\x8c\t\x8d\xf8\x10\xad\x8aD\xcc\x9cc\xb5ZR\xbc>6\x169\xect\x99\xa8Dn\x82\xff\x99\x1f\x07\xccX\xdf\xc0\xe5\x9e\x8dR\xc7\x89\xa98\xf2S\x86\xddI\xa6\x89\xd3rO\xf1\xeb\xf7&\xe8*\x93L^#\xa5\x95p\x99\xcb\xb0\xb2\xad\xf0\xfc\x80\xd2\xe2\xcd\x01\xb7\xccF\xda\x04\x8ed\x19K\xb4\x17\x81\x1aZ*/k\x0c\xdb\xb2\x03[\x14\x03-j\x0c%\xe3\xe6\xa3\x0e\xb61\xcf\xbe\x8f\x18\x06\x85\xf8\xf4r\xc0r \x03j\x02\xe9\x0f\xc4r\xe1\x01xw\xbe>6~\x84\xe9\x8d~N\xdc/\xfb\xb2\x18C\xfd\x9b\xfe\x81~\xa2\xfe\xef\xcew \xb8\x02\xe4\xb9KV]\xb7\xbb\xe0\r\xe1\x8f|\x84\xc0a\x1bZ8\xb3E\xc4\xd0M\xed\xa4\xd7\x92gR.f\xae\xb6\xc1`rT\xd6\xd4\x96-\x05\x0e3JL\xcf\xfe*\xf7lAx\x85Fgnz\x8byc^]\x9a\r\xce\xab\xc1/\xab\x0cG\xd0\xa9\xec\x94\xb2\xe6\x16\xfa\xfc\xc1M\xd2Xk\x0c>\xc4\x0f\xcd+\x82\xb7\xdbs5\x03\xf8\xc8vF\xb9FI\xfe\xc7T|\xbb\xfe1@I\x12\\#h\x96\x04\xf2\xb0^\xd4\x8c\xdaL\x1b\x1a\x83\x0f\xf9i\xca\xc5\xea\x01\x1f<,\xbc;6`R\xebf`Y]\xff\xc5<?>\xc9\xed\xca\xaf\'\x9d"\x83\x03/G!|\x86\xd8\xb7\x15P\x80\x14\x0e(?\xfe\x1c\xf9\x92f5\x1f\xb9\xa4\x17\xad\xba\xa9Ev\x84pr\x85\x99\x18\x97\x14\xaa\x8dg\x05\xaa\x0cop3\xa7\xa53\xe7\x1f5\x8a\x86Y\xe5:&\xad\x05X\xee=ev\x93\xb3\x06\xb0RJA\x17Rz\x12\xd8u\x93n}\x15P\xca\xc9\x08\xdb/\xd4\xecz\xbc\xf5\xba\xf7a\x9b\x81H\xda\x1aV\x08\xa2\n\x18\xcf\xe1/\x87L\xf6_7\x8c\xed\xfc\xeez\x16\xb2y\xe1\xc9(D\x175E\x91\xb2\xfe\xbb\xb1l=\xa3\x0f5s`\x1d\x0f\x04\x96&\xe0\x91@7\x91\x94\xf8\x1d\xf4\xccX\xd3k\xf4f\xfb\x81\xfb\n\r\xc8\x02&#n\x82)\xa7\xc4x\x8c\x161|z3\x0f\xeb\xd3g<\x86!f \x13K\xab\xcb\x93zqM\x03]\xc9PDaWbe\x1c\x01\xa6\x97\n \xd4\x8eKD\xbf]S\x81\x8f\xd3y\xd2\x08Z\xd9\xd9wc\x93\xb2]\xa7\x1c\x94\xe1e\xf0\x92\x9e\xe1\xb3r\xfe\x80\xf3J\xad\x0beq\x94\xfe\xa3\x15\x85%Z\x00\xe1\xc227\x95/\xfe\x0b\xf2\xfe\x93xE\'"\x9f\xca\xe1\x17\x8ds\xfe\xdc\xc1z\xe7Q\xca\xee\xb3\xf0\x845\xde3\x87C1\xc9\x0e\xa6\xfc\xf7\xa4\xa1#\x9fg\x12\xf5\xe7\xffG\xd7\x8b\xc1n6\x9aQ\xbc\x9fhY\xe6C\xd2\x16\x8b\xa0\x84)\xb8Y\x1b\x1d\xf0x\xd9*\x1f\xe3\x16\x86)\xa4\x91\xbfS\xea\xecV\xef\xec\x04\xc8\xa7\xc6g\x9eq\xe9m\xc4\xc0\r\x81f\x86}GNwq_KzX\xab\x16\xf7\xf6\xfbS\xe0\x97\xc9\x95/1a\x12\x01k\x0f\xe7\xed\xe4\xe6\x15\xb8p\xd1\xbd\xe9\xc1Nn\x96\x84y[\xb54?\x00\x86N\x1e\xeb.\xea\xd0\xae\x068\xc2~\x1fj6\xe5i\xb3\x9a\xbc6\xb6\x95\xcc\xf7>GZ\xf9v\x8c\x9f\x17_@\x9e8\xcd\x19\xd9\xa0\xe0\xcc\xb1\xf3\x05\x90\xeaKu\x19>\x9cYN#)\x0b\xa7\x97\xae?\x13\x0f\xe7\xd3\x7f!\x1b\xecG\xa7\xf2\'\x90\r\x87\xee\xee[\xd3vF\x88\xca1\xb8\xc9\x1d\xd3p\xb1\xe7\xad\xd09\x06\x1f\x02\x00W\xd9\x9b\x99d"fnz\x90;\xe3\x03\x06\x9b\x92\xc3\x03xq\x97g\xe6\xadR\xce(\xf9\x8c\xa4N\x86\x1fB+Y\xcd\xc1\xba\xbc\xa2\x1c8\xd2\xb1\xccP\xe2\\\xfa\x9e@\xca\xfbS\x9fp\xba\xab\xa8\xee\xe0\x8e\xdf\xbc#\xacd\x1aP\x1c4\xady\x17N\xe3\xc9\xe9<\xf0\xcfR\xb6\xa6p\xcf\xc8\x16[v\x9d"\x1fS\xebhCZ\x1f3\x13_\xc8\xb5\xb0\xaf\xe7F\xd9\xb98\x1a\x0c:`P\x01\x8b ^\xca\x94w7\xd5\xb5\xd2\x88|K\r\xa7:\x00\xe0\t\xedy\xa1\xc80\xd0g{\xa5\x85\xed33%J9\xdf\x19\x11rk\x0b]3Q\xda\xd0\x14\xe9\x98\xaeD&\x1b\xb9\x00\xb8\x9b-gI=\x19\xbf4E,\xfb<\x10\xf5\xf66/\x08\xcdy~&9\xfe\x8cf\xfeg\x13\x1dM]\x94\xcc\x17\xadid\xd6\xcb\xfaU\xf0\x95\xa5\xa1\x19<j\xaf\xe1z\xdb\xfe\x94\x9a\x85y\xcc\x14^\x8d\x1b\xd5\x0c\x1a\xe2z/\xf9\xe1\t/\xa4I\x80w\x8c\xcd\xe0\x97\xddi\x90\xed\xe0\x04-\xa3\x06\xc2\xddv9\xff\xb2l\x1b\\G\xd7cL1\xcb\x07\xafm\xb7)\xf7\x08\xef\xe7\xac7\xcd\x04\x0cd\x0cP\xff\x06\x97WFv\x8c\x9e\x054\xe09\xbd\xe5Io%\xf54\x9d\xab5\xea\xd1\xcb\xfe\xbe\xe11<B\xb6\x03y5\xe6\x98\xda<\x17gR\x84\x94J\xd3*E)\xab\xffG\xd5\xf3\xab\xaf\x02\xf9hE+\xc8\xd5\x87\xad\x95.N\xf1~\xca\xda\xd8MJ\x94\x1a\x02i\xe3\x8c\x11\xb0oU\xa6\x04\xb0\xdb\xd0D\'\xe7\x83\xb6\xd8\xad\xdd\x8d\xe65 u\xd3\xb8\x0b\xcc\xe1\xf1\x8av\x1a\x05\xee\xf1\x02\xbfz\xed1\x8a\xbb\xa1)\xe2P<\x81z\x14\x92/\n\xc9\xc8\x82\x9c\x85\xc4tc\x82[\xa2Q\xbb\x99\xc0\x9aP\xddZ\x11\xde\x0c\xb4\xd0\xcbl\xf6\x1e\x9b\x90S&}%m(\xed@\xcff\xea1\'\x83H\x93\xcb|f\xad\x8bvt\x80\xe5\x84\xac\x92\xee\x8d$\xc3\x8b\xb8\xc5X\xf5\x96O\xe1\x83L\x8b\xf3\x00\xb3\xf2\x13\xa2*J\xeb\x15\x85\x19*\x7f\xd3\xf1w$\xee*\xaf=<\x0cb\xe4c\x16\xa96\xe691\x06\xcd\xcc4\x0e\x9f_\xc4\x18\xb9\x1aEm\xcc\x01\x98\xb4Ns2\xc7WAg\x04\x8eM*\xc5\xc9\x8d\xbb\xfd\x84q&\x94V\xa0\xb4\xd74\x8ea\x9dY\x8bo\x86\xf6\xdeH\xcd\xff\xf4\xe5\xd9\xe5H\xc8ML\xf8\x94\x80\xb8\x88\x00\x82\xee\xe7\x9eN\x1f\xaeA\x85\xd9\x8f\x8b\x87\x8d\x94Tu\n.\xc0\x9d^\xe4\x8e\nt\xbbP\x93\xea\x01\xd9C|L\xdc\xf4\xcd\xaf\xbet\x95\x90\xafy\xff\x9b6\x91\x87\x1d\xa9\x16l2\x8d\xf8\xf0\x07Z\xb5\x938]\x0b\x1e\r\x1d\xd9\xea\x06G\xae\x06\xa0'
    tDesObj = TDES(secret_key3, "OFB", initvector3)
    print(tDesObj.decrypt(ciphertext3).decode("utf-8"))