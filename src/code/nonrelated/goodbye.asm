section .data
    helloworld db "Hello, World!", 0x0A    

section .text
    global main

main:


    mov rax, 1                  ; syscall number for sys_write (1)
    mov rdi, 1                  ; file descriptor 1 is stdout
    lea rsi, [rel helloworld]   ; load effective address of helloworld into rsi
    mov rdx, 14                 ; hardcoded number of bytes to write (14)
    syscall                     ; invoke the syscall to write to stdout



    mov rax, 60                 ; syscall number for sys_exit (60)
    xor rdi, rdi                ; set exit code to 0 (success)
    syscall                     ; invoke the syscall to exit the program
