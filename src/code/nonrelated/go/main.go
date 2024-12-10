package main

import (
	"fmt"
	"sync"
)

func display(str string) {
	for i := 0; i < 3; i++ {
		fmt.Println(str)
	}
}

func main() {
	var wg sync.WaitGroup

	wg.Add(1)
	niga := func() {
		display("Hello, niga!")
		wg.Done()
	}
	go niga()

	wg.Add(1)
	go func() {
		display("Goodbye, nige!")
		wg.Done()
	}()

	wg.Wait()

	display("Hello, Main!")
}
