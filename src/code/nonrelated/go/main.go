package main

import (
	"fmt"
	"strings"
)

func main() {
	var t int
	fmt.Scan(&t)

	for i := 0; i < t; i++ {
		var str string
		fmt.Scan(&str)

		new := strings.Repeat(" ", len(str))
		index := 0

		for i := len(str); i > 0; i-- {
			switch(str[i]) {
			case: 'w'
			break;

			case: 'q'
				new[index++] = 'p'
			break;

			case: 'p'
				new[index++] = 'q'
			break;

			default:
				break;
			}
		}

		fmt.Println(new)
	}
}
