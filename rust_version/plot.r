# Числовые параметры задачи
source("ESC_parameters.r")
# Параметры расположения наноструктур (из файла)
nanostructures <- read.table("results/nanostructures.dat")
xn <- nanostructures$V1
yn <- nanostructures$V2
cn <- nanostructures$V3
# Матрица с номерами ячеек
cells <- matrix(1:Nc, Ncx, Ncx)
# Функция для определения номера ячейки, куда попадает данная точка
cell_numbers <- function(x, y){
	cx <- floor(x / dc) + 1
	cy <- floor(y / dc) + 1
	if(cx > Ncx){cx <- Ncx}
	if(cy > Ncx){cy <- Ncx}
	cyd <- if(cy == 1){Ncx}else{cy - 1}
	cyu <- if(cy == Ncx){1}else{cy + 1}
	cxl <- if(cx == 1){Ncx}else{cx - 1}
	cxr <- if(cx == Ncx){1}else{cx + 1}
	c(cells[cx, cy], cells[cx, cyu], cells[cx, cyd], cells[cxl, cy], cells[cxr, cy], cells[cxl, cyu], cells[cxl, cyd], cells[cxr, cyu], cells[cxr, cyd])
}

# Посчитаем, сколько центров наноструктур находится в каждой ячейке
nano_in_cell <- rep(0, Nc)
# Создадим список с координатами
nano_list <- matrix(rep(0), Nc, N)
# Вот теперь распределяем наноструктуры по ячейкам
for(jn in 1:N){
	jc <- cn[jn]
	nano_in_cell[jc] <- nano_in_cell[jc] + 1
	nano_list[jc, nano_in_cell[jc]] <- xn[jn] + 1i * yn[jn]
}
#===============================
# Функция, задающая вид наноструктуры
U1 <- function(x, y, xn, yn){
        if(abs(x - xn) < 5 * a && abs(y - yn) < 5 * b){
            exp(-(x - xn)^2 / a^2 / 2 - (y - yn)^2 / b^2 / 2)
        } else{0}
    }
# Общая функция для потенциала всех наноструктур
U <- function(x, y){
	nc <- cell_numbers(x, y)
	Ures <- 0
	for(jc in nc){
    sc <- nano_in_cell[jc]
	if(sc > 0){
		for(j in 1:sc){
			z <- nano_list[jc, j]
			Ures <- Ures + U1(x, y, Re(z), Im(z))
		}
		}
	}
    Ures
}
# Параметры сетки для построения графика потенциала (это нужно только для иллюстрации)
dx <- dy <- 3
Nx <- Ny <- floor(L / dx)
Ugrid <- matrix(rep(0), Nx, Ny)
# Построение графика
for(jx in 0:(Nx-1)){
    for(jy in 0:(Ny-1)){
        x <- jx * dx
        y <- jy * dy
        Ugrid[jx+1,jy+1] <- U(x, y)
    }
}
# Считываем данные для 5 траекторий
trajectories <- matrix(rep(0), 10, Nt - 1)
velocities <- matrix(rep(0), 5, Nt - 1)
for(je in 1:5){
    trajectory <- read.table(paste("results/Trajectory_", je, ".dat", sep = ""))
    trajectories[2*je-1,] <- trajectory$V2
    trajectories[2*je,] <- trajectory$V3
    velocities[je,] <- trajectory$V4
}
# Экспорт изображения наноструктур с траекториями нескольких электронов
flnm <- "Trajectories.png"
png(flnm, width = 3000, height = 3000)
par(mar = c(18, 22, 5, 5))
par(mgp = c(12, 4, 0))
image(1:Nx * dx, 1:Ny * dy, Ugrid, col = hcl.colors(100, "RdBu", rev = TRUE),
xlab = "x, nm", ylab = "y, nm", cex.axis = 6, cex.lab = 8)
points(trajectories[1,], trajectories[2,], col = "black", pch = 16, cex = 2)
points(trajectories[3,], trajectories[4,], col = "green", pch = 16, cex = 2)
points(trajectories[5,], trajectories[6,], col = "purple", pch = 16, cex = 2)
points(trajectories[7,], trajectories[8,], col = "red", pch = 16, cex = 2)
points(trajectories[9,], trajectories[10,], col = "orange", pch = 16, cex = 2)
dev.off()
#=========================
# Теперь проанализируем среднюю скорость за все время движения для всех электронов
vx_mean <- rep(0, Ne)
for(je in 1:Ne){
    trajectory <- read.table(paste("results/Trajectory_", je, ".dat", sep = ""))
    vx_mean[je] <- mean(trajectory$V4)
}
# Построим гистограмму
flnm <- paste("Velocities_Mean.png");
png(flnm, width=3000, height=2500)
par(mar = c(23,22,5,5))
par(mgp = c(12,4,0))
hist(vx_mean, xlab = "Средняя скорость, нм/фс", ylab = "Количество электронов", cex.lab=8, cex.axis=7, col = "skyblue3", breaks = 12)
box(lwd=4)
dev.off()
# Выводим все средние скорости в один файл
write.table(vx_mean, "Velocities_Mean.dat", col.names = FALSE, row.names = FALSE)
vx <- read.table("Velocities_Mean.dat")
paste("Средняя скорость равна ", mean(vx$V1), "нм/фс")
paste("при напряженности поля ", eE, "В/нм")
# Конец программы
