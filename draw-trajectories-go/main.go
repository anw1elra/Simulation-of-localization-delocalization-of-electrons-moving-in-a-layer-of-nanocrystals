package main

import (
	"bufio"
	"fmt"
	"github.com/fogleman/gg"
	"os"
)

func main() {
	im, err1 := gg.LoadImage("potential_map.png")
	if err1 != nil {
		panic("Error opening potential_map.png file")
	}

	var nE int

	//fmt.Print("Enter the largest E index: ")

	//fmt.Scan(&nE)

	nE = 50

	for jE := 1; jE <= nE; jE++ {
		dc := gg.NewContextForImage(im)
	
		for je := 1; je <= 10; je++ {
			filepath := fmt.Sprintf("results/Trajectory_e%v_E-%v.dat", je, jE)
			file, err2 := os.Open(filepath)
			if err2 != nil {
				panic("Can't open trajectory file")
			}
			defer file.Close()

			scanner := bufio.NewScanner(file)

			for scanner.Scan() {
				var x, y float64
				fmt.Sscanf(scanner.Text(), "%f %f", &x, &y)
				
				switch je {
				case 1:
					dc.SetRGB255(31,119,180)
				case 2:
					dc.SetRGB255(255,127,15)
				case 3:
					dc.SetRGB255(44,160,44)
				case 4:
					dc.SetRGB255(214,39,40)
				case 5:
					dc.SetRGB255(148,103,189)
				case 6:
					dc.SetRGB255(140,86,75)
				case 7:
					dc.SetRGB255(227,119,194)
				case 8:
					dc.SetRGB255(118,183,178)
				case 9:
					dc.SetRGB255(237,201,73)
				case 10:
					dc.SetRGB255(23,190,206)
				}

				for jx := -2; jx <= 2; jx++ {
					for jy := -2; jy <= 2; jy++ {
						dc.SetPixel(int(x)+jx, int(y)+jy)
					}
				}
			}
			
		}
		dc.SavePNG(fmt.Sprintf("plots/Trajectory_E-%v.png", jE))
	}
}
