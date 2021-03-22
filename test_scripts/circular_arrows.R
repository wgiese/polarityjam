library(ggplot2)

bins <- 10
bin_size = 360/bins

angle_deg <- c(0,20,30,20,100,180,23)
pdf(paste0(getwd(),"test.pdf"))

p <- ggplot() +
	geom_histogram(aes(x=angle_deg, y = ..ncount..),
	     breaks = seq(0, 360, bin_size),
	     colour = "black",
	     fill = "grey80") +
	ggtitle("cellular orientation") +
	theme(axis.text.x = element_text(size = 18)) +
	coord_polar(start = -pi/2.0, direction = -1) +
	scale_x_continuous(limits = c(0, 360),
		 breaks = (c(0, 90, 180, 270))) +
	theme_minimal(base_size = 14) +
	ylab("") +
	theme(axis.text.y=element_blank()) +
	scale_y_sqrt()
p <- p + geom_segment(aes(x=100, y=0, xend=100, yend=1, color="red"))
p

dev.off()
