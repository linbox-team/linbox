// Anti-Spam Script Adapted from Laurent Fousse (c) 2007
// Time-stamp: <11 Oct 11 18:41:02 Jean-Guillaume.Dumas@imag.fr>
function CryptEMail(Address) {
    cipher = "<f mwjk=\"rfnqyt:" + Address + "\">" + Address + "</f>"
    alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key = 5
    clear = ""
    for (i = 0; i < cipher.length; i++) {
	if (alpha.indexOf (cipher.charAt (i)) == -1) {
	    clear += cipher.charAt (i)
	} else {     
            if (alpha.indexOf (cipher.charAt (i)) < key) {
	    offset = (alpha.indexOf (cipher.charAt (i)) + 26 - key) % alpha.length
	    clear += alpha.charAt (offset)
            } else 
            {
            
	    offset = (alpha.indexOf (cipher.charAt (i)) - key) % alpha.length
	    clear += alpha.charAt (offset)
            }
            
	}                               
    }
    document.write (clear)
}

function CryptedMail(Address,Name) {
    cipher = "<f mwjk=\"rfnqyt:" + Address + "\">"
    alpha = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key = 5
    clear = ""
    for (i = 0; i < cipher.length; i++) {
	if (alpha.lastIndexOf (cipher.charAt (i)) == -1) {
	    clear += cipher.charAt (i)
	} else {     
	    offset = (alpha.lastIndexOf (cipher.charAt (i)) + alpha.length - key) % alpha.length
	    clear += alpha.charAt (offset)
	}                               
    }
    document.write (clear + Name + "</a>")
}

function emailcrypted(Address) {
    cipher = "rfnqyt:" + Address
    alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key = 5
    clear = ""
    for (i = 0; i < cipher.length; i++) {
	if (alpha.indexOf (cipher.charAt (i)) == -1) {
	    clear += cipher.charAt (i)
	} else {     
	    offset = (alpha.indexOf (cipher.charAt (i)) + alpha.length - key) % alpha.length
	    clear += alpha.charAt (offset)
	}                               
    }
    window.location.href = clear;
}
