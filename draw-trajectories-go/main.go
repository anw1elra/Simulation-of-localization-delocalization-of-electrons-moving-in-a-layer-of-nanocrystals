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
					dc.SetRGB255(255,56,24)
				case 2:
					dc.SetRGB255(129,68,163)
				case 3:
					dc.SetRGB255(255,125,36)
				case 4:
					dc.SetRGB255(255,217,57)
				case 5:
					dc.SetRGB255(163,81,58)
				case 6:
					dc.SetRGB255(175,211,80)
				case 7:
					dc.SetRGB255(117,213,209)
				case 8:
					dc.SetRGB255(78,190,98)
				case 9:
					dc.SetRGB255(56,204,246)
				case 10:
					dc.SetRGB255(16,123,196)
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
