import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import os

def load_density_file(filename, grid_size):
    # Read the raw binary file
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
        data = data[3:] # Skip the first 3 elements
    
    # Reshape the 1D array into a 3D grid with x-major ordering
    return data.reshape((grid_size, grid_size, grid_size), order='F')  # Fortran-style (column-major) ordering

def visualize_density_slices(filename, grid_size=128, output_dir='plots'):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the density data
    density = load_density_file(filename, grid_size)
    
    # For small grids (3x3 or smaller), show all slices
    if grid_size <= 3:
        # Calculate number of slices in each dimension
        n_slices = grid_size
        
        # Create a figure with ImageGrid for all slices
        fig = plt.figure(figsize=(15, 5 * n_slices))
        grid = ImageGrid(fig, 111,
                        nrows_ncols=(n_slices, 3),
                        axes_pad=0.1,
                        share_all=True)
        
        # Plot all slices
        for i in range(n_slices):
            # XY slice (top view)
            im0 = grid[i*3].imshow(density[:, :, i], cmap='viridis')
            grid[i*3].set_title(f'XY Slice (z={i})')
            
            # XZ slice (front view)
            im1 = grid[i*3+1].imshow(density[:, i, :], cmap='viridis')
            grid[i*3+1].set_title(f'XZ Slice (y={i})')
            
            # YZ slice (side view)
            im2 = grid[i*3+2].imshow(density[i, :, :], cmap='viridis')
            grid[i*3+2].set_title(f'YZ Slice (x={i})')
        
        # Add colorbar
        plt.colorbar(im0, ax=grid.axes_all)
        
    else:
        # For larger grids, show only center slices
        fig = plt.figure(figsize=(15, 5))
        grid = ImageGrid(fig, 111,
                        nrows_ncols=(1, 3),
                        axes_pad=0.1,
                        share_all=True)
        
        # Plot three orthogonal slices through the center
        mid = grid_size // 2
        
        # XY slice (top view)
        im0 = grid[0].imshow(density[:, :, mid], cmap='viridis')
        grid[0].set_title('XY Slice (Top View)')
        
        # XZ slice (front view)
        im1 = grid[1].imshow(density[:, mid, :], cmap='viridis')
        grid[1].set_title('XZ Slice (Front View)')
        
        # YZ slice (side view)
        im2 = grid[2].imshow(density[mid, :, :], cmap='viridis')
        grid[2].set_title('YZ Slice (Side View)')
        
        # Add colorbar
        plt.colorbar(im0, ax=grid.axes_all)
    
    # Get frame number from filename
    frame_num = filename.split('_')[-1].split('.')[0]
    
    plt.suptitle(f'Density Visualization - Frame {frame_num}')
    
    # Save the plot
    output_filename = os.path.join(output_dir, f'slices_{frame_num}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved slices visualization to {output_filename}")

def visualize_volume(density_data, step, output_dir='plots'):
    """Create a 3D volume visualization of the density field."""
    # Create figure with 3D projection
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create a 3D grid - note the order of dimensions
    x, y, z = np.meshgrid(
        np.arange(density_data.shape[0]),
        np.arange(density_data.shape[1]),
        np.arange(density_data.shape[2]),
        indexing='ij'  # Use 'ij' indexing to match the data layout
    )
    
    # Plot voxels where density is above threshold
    threshold = 0
    mask = density_data > threshold
    
    # Plot the voxels - ensure mask is boolean
    ax.voxels(mask, facecolors=plt.cm.viridis(density_data), alpha=0.3)
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Fluid Simulation - Step {step}')
    
    # Set the viewing angle
    ax.view_init(elev=30, azim=45)
    
    # Save the figure
    filepath = os.path.join(output_dir, f'voxels_{step:04d}.png')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close(fig)

def process_all_frames(start_frame=0, end_frame=100, grid_size=128):
    """Process all frames in the range and save visualizations"""
    for frame in range(start_frame, end_frame):
        filename = f"./renderData/density_{frame:04d}.bin"
        velocity_filename = f"./renderData/velocityX_{frame:04d}.bin"
        if os.path.exists(filename):
            print(f"Processing frame {frame}...")
            visualize_density_slices(filename, grid_size)
            #visualize_density_slices(velocity_filename, grid_size)
            #visualize_volume(load_density_file(filename, grid_size), frame)
        else:
            print(f"Warning: File {filename} not found")

if __name__ == "__main__":
    # Create plots directory
    os.makedirs('plots', exist_ok=True)
    
    # Process all frames from 0 to 99
    process_all_frames(0, 100, 64)
    
    print("\nAll visualizations have been saved to the 'plots' directory")
    print("Files are named as:")
    print("- slices_XXXX.png for 2D slice visualizations")
    print("- volume_XXXX.png for 3D volume visualizations")